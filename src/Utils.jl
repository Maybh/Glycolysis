using Distributed, ClusterManagers, Dates, CMAEvolutionStrategy, JSON, CSV

function save_csv_with_subdir(data, dir, filename)
    if !isdir(dir)
        mkdir(dir)
    end
    # Define the file path within the subdirectory
    filepath = joinpath(dir, filename)

    # Write the DataFrame to the CSV file
    CSV.write(filepath, data)
end

"load configuration from json file"
function read_config(path)
    config = open(path, "r") do file
        JSON.parse(file)
    end

    config["param_names"] = Symbol.(config["param_names"])
    config["metab_names"] = Symbol.(config["metab_names"])
    config["nt_param_choice_names"] = Symbol.(config["nt_param_choice_names"])
    config["non_opt_params"] = Symbol.(config["non_opt_params"])

    for (name, value) in pairs(config["rescaling"]) 
        config["rescaling"][name] = eval(Meta.parse(value)) 
    end

    return config
end
 
"data loader with preparation"
function load_prepare_data(file_path)
    #Load enzyme kinetic data data
    data = CSV.read(file_path, DataFrame)

    #Add source column that uniquely identifies a figure from publication
    data.source .= data.Article .* "_" .* data.Fig

    data.fig_num = vcat(
        [
            i * ones(Int64, count(==(unique(data.source)[i]), data.source)) for
            i = 1:length(unique(data.source))
        ]...,
    )
    
    return data
end

"rescaling parameters by given function, where {x} will replace with the respective value"
function param_rescaling_from_conf(nt_params_names, rescaling_conf=Dict())
    nt_params_names_new = Dict{Symbol, Float64}()
    for (name, value) in pairs(nt_params_names)
        if !isempty(rescaling_conf) && string(name) in keys(rescaling_conf)
            nt_params_names_new[name] = rescaling_conf[string(name)](value)
        else
            nt_params_names_new[name] = value
        end
    end
    
    return NamedTuple{Tuple(keys(nt_params_names))}((nt_params_names_new[n] for n in keys(nt_params_names)))
end

function param_subset_select_from_conf(nt_params_names, nt_param_choice, param_constraints, param_names, constant_params)
    nt_params_names_new = Dict(pairs(nt_params_names))
    
    # set costant value for parameters
    for (name, value) in pairs(constant_params)
        nt_params_names_new[Symbol(name)] = value
    end
    
    # set value for K parameters, 1 -> K_i=K_a, 2 -> K_a=Inf and 3-> K_i=Inf
    for (name, choice) in pairs(nt_param_choice)
        name_str = string(name)
        choice_str = string(choice)

        if startswith(uppercase(name_str), "K")
            if choice > 0
                for (param, value) in pairs(param_constraints[name_str][choice_str])
                    if choice == 1
                        nt_params_names_new[Symbol(param)] = nt_params_names_new[Symbol(value)] 

                    elseif choice == 2
                        nt_params_names_new[Symbol(param)] = Inf
        
                    elseif choice == 3
                        nt_params_names_new[Symbol(param)] = Inf
                    end
                end
            end

        # Then, check if it's a dictionary and the choice exists
        elseif haskey(param_constraints, name_str) && isa(param_constraints[name_str], Dict) && haskey(param_constraints[name_str], choice_str)
            choice_dict = param_constraints[name_str][choice_str]

            # Iterate over the inner dictionary
            for (param, value) in pairs(choice_dict)

                # Check the type of value and handle accordingly
                if isa(value, Float64) || isa(value, Int)
                    nt_params_names_new[Symbol(param)] = value

                elseif isa(value, String)

                    # Handle string values (e.g., "Inf" or parameter names)
                    if value == "Inf" || value == "inf"
                        nt_params_names_new[Symbol(param)] = Inf

                    else
                        nt_params_names_new[Symbol(param)] = nt_params_names_new[Symbol(value)]
                    end
                end
            end       
        end
    end

    return NamedTuple{Tuple(param_names)}((nt_params_names_new[n] for n in param_names))
end


"initial slurm cluster"
function initial_cluster()
    #Increase time before worker terminates without hearing from master process. Useful for large number of cores.
    ENV["JULIA_WORKER_TIMEOUT"] = 600.0

    # ENV["SLURM_JOB_CPUS_PER_NODE"] give number of cores in "40(x4),32,20" format
    # "40(x2),32,20" can be parsed to 132 using code below to get total number of allocated cores

    subs = Dict("x" => "*", "(" => "", ")" => "");
    np = sum(eval(Meta.parse(replace(ENV["SLURM_JOB_CPUS_PER_NODE"], r"x|\(|\)" => s -> subs[s]))))
    addprocs(SlurmManager(np); exeflags = "--project")
end

function initial_local_distribution(;n_process=8)
    #Uncomment code below to run on local machine or uncomment code above to run on cluster
    #If run locally adjust number of workers to number of cores on your machine
    addprocs(n_process; exeflags = "--project")
end

function filter_out_subsets_nan_values(data, param_subsets, config)
    param_names = config["param_names"]
    metab_names = config["metab_names"]
    nt_param_choice_names = config["nt_param_choice_names"]

    # convert DF to NamedTuple for better type stability / speed
    rate_data = Tables.columntable(data[.!isnan.(data.Rate), [:Rate, metab_names..., :fig_num]],)
    # make a vector containing indexes of points corresponding to each figure
    fig_point_indexes = [findall(data.fig_num .== i) for i in unique(data.fig_num)]


    param_subsets = [
        NamedTuple{Tuple(nt_param_choice_names)}(x) for
        x in unique(param_subsets) ]

    # Filter out parameter subsets that produce NaN loss values (e.g., when all substrate K = Inf)
    param_subsets = [
        param_subset for param_subset in unique(param_subsets) if !isnan(
            loss(
                ones(length(param_names)),
                rate_data,
                fig_point_indexes, 
                config;
                nt_param_choice = param_subset,
            ),
        )
    ]
    return param_subsets
end

function redirect_to_files(dofunc, outfile, errfile)
    open(outfile, "w") do out
        open(errfile, "w") do err
            redirect_stdout(out) do 
                redirect_stderr(err) do 
                    dofunc()
                end
            end
        end
    end
end
