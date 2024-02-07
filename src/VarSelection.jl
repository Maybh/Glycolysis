using Distributed, DataFrames
include("Utils.jl")
include("FittingFunction.jl")
include("Enzyme.jl")

function create_param_subsets(param_constraints, nt_param_choice_names)
    # Make an array of all possible param combination
    combination = Vector{Any}()
    for name in nt_param_choice_names
        name_str = string(name)
        if haskey(param_constraints, name_str) 
            constrains = param_constraints[name_str]
            n_choices = haskey(constrains, "0") ? length(constrains) - 1 : length(constrains)

            push!(combination, 0:n_choices)
        else
            push!(combination, 0)
        end
    end 

    return collect(Iterators.product(combination...))
end

function find_best_models_per_fig_all_subsets(
    data,
    all_subset_filterd, 
    config;
    n_non_params=2, 
    n_iter = 20, 
    maxiter = 50_000, 
    save_results = false
    )

    figs = unique(data.source)
    all_subsets_figs_to_fit = collect(Iterators.product(all_subset_filterd, figs));
    
    results_array = pmap(
        subset_fig_to_fit -> fig_loocv(subset_fig_to_fit[2], 
        data, 
        config; 
        n_iter = n_iter, 
        nt_param_choice = subset_fig_to_fit[1], 
        maxiter = maxiter),
        all_subsets_figs_to_fit,
    )

    results_df = DataFrame(
        n_removed_params = Int[], 
        error_counter = Int[],
        param_subset = String[], 
        dropped_fig = String[], 
        params = Vector{Float64}(undef,0), 
        train_loss = Float64[], 
        test_loss = Float64[]   
    )

    for i in eachindex(results_array)
        param_subset = all_subsets_figs_to_fit[i][1]
        keys_to_exclude = config["non_opt_params"]
        n_removed_params = sum((param_subset[key] > 0) for key in fieldnames(typeof(param_subset)) if !(key in keys_to_exclude))
        params = NamedTuple{Tuple(config["param_names"])}(results_array[i][4])

        row_dict =  Dict("n_removed_params" => n_removed_params,
            "error_counter" => results_array[i][5],
            "param_subset" => string(param_subset),
            "dropped_fig" => results_array[i][1],
            "train_loss" => results_array[i][2],
            "test_loss" => results_array[i][3],
            "params" => params,  
        )

        results_df = push!(results_df, row_dict, promote=true)  
    end
    
    # Filter subsets with NaN train or test loss values
    filter!(row -> !isnan(row.train_loss), results_df)
    filter!(row -> !isnan(row.test_loss), results_df)
    
 #     # these lines creates a column for each param: 
 #     max_length = maximum(length.(df_results_all.params))
 #     results_params_df = DataFrame([getindex.(df_results_all.params, i) for i in 1:max_length], :auto)
 #     rename!(results_params_df, param_names)
 #     df_results_all = hcat(select(df_results_all, Not(:params)), results_params_df)   
    
    # keep for each fig and n removed param: only subsets with training loss <= 1.1 * minimal training loss
    # Group by 'dropped_fig' and 'n_removed_params' and calculate the minimum 'train_loss' for each group
    min_train_loss = combine(groupby(results_df, [:dropped_fig, :n_removed_params]), :train_loss => minimum => :min_train_loss)
    # Merge the minimum train loss back into the original data frame
    results_df = leftjoin(results_df, min_train_loss, on = [:dropped_fig, :n_removed_params])
    # Filter rows where 'train_loss' is less than 1.1 times the minimum for the group
    df_best_subsets = results_df[results_df.train_loss .< 1.1 .* results_df.min_train_loss, :]

    # Drop the 'min_train_loss' column
    select!(df_best_subsets, Not(:min_train_loss))
    
    return df_best_subsets
end


function find_best_subsets_per_fig(dropped_fig, all_param_subsets,data, config; 
    n_iter=20, n_non_params=2, maxiter = 50_000, print_prog = print_prog)
    
    param_names = config["param_names"]
    metab_names = config["metab_names"]
    nt_param_choice_names = config["nt_param_choice_names"]

    if print_prog ==true
        println("figure: ", dropped_fig)
    end
    # convert DF to NamedTuple for better type stability / speed
    rate_data = Tables.columntable(data[.!isnan.(data.Rate), [:Rate, metab_names..., :fig_num]],)
    
    # make a vector containing indexes of points corresponding to each figure
    fig_point_indexes = [findall(data.fig_num .== i) for i in unique(data.fig_num)]
    
    df_results_all = DataFrame(
        n_removed_params = Int[], 
        error_counter = Int[],
        param_subset = String[], 
        dropped_fig = String[], 
        params = Vector{Float64}(undef,0),
        train_loss = Float64[], 
        test_loss = Float64[]   
    )

    global starting_param_subsets = all_param_subsets[0]
    # Initialize array to save the best models from each round
    best_param_subsets = NamedTuple[]

    for n_removed_params = 0:7
        starting_param_subset_masks = unique([
            (
                mask = [(starting_param_subset[1:end-n_non_params] .== 0)..., zeros(Int64, n_non_params)...],
                non_zero_params = starting_param_subset .* (starting_param_subset .!= 0),
            ) for starting_param_subset in starting_param_subsets
        ])
        
        subset_of_all_param_subsets_w_n_removed_params = all_param_subsets[n_removed_params]
   
        
        # actually choose param_subsets with n_removed_params number of parameters removed that also contain non-zero elements from starting_param_subsets
        param_subsets = []
        for starting_param_subset_mask in starting_param_subset_masks
            push!(
                param_subsets,
                unique([
                    subset_w_n_removed_params .* starting_param_subset_mask.mask .+
                    starting_param_subset_mask.non_zero_params for
                    subset_w_n_removed_params in subset_of_all_param_subsets_w_n_removed_params if sum(
                        (
                            subset_w_n_removed_params .* starting_param_subset_mask.mask .+
                            starting_param_subset_mask.non_zero_params
                        )[1:end-n_non_params] .> 0,
                    ) == n_removed_params
                ])...,
            )
        end
        
    
        param_subsets = [
            NamedTuple{Tuple(nt_param_choice_names)}(x) for
            x in unique(param_subsets) if sum(values(x[1:end-n_non_params]) .> 0) == n_removed_params
        ]
    
                
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
        @info "n_removed params: $n_removed_params, n param_subsets: $(length(param_subsets))"
        @info "starting fig_loocv, $dropped_fig"
        results_array = pmap(
            param_subset -> fig_loocv(dropped_fig, data, config; maxiter = maxiter, n_iter = n_iter, nt_param_choice = param_subset),
            param_subsets,
        )
        @info "finished fig_loocv, $dropped_fig"        
        results_df = DataFrame(
            n_removed_params = Int[], 
            error_counter = Int[],
            param_subset = String[], 
            dropped_fig = String[], 
            params = Vector{Float64}(undef,0), 
            train_loss = Float64[], 
            test_loss = Float64[]   
        )

        for i in eachindex(results_array)
            # params = NamedTuple{Tuple(param_names)}(results_array[i][4])

            row_dict =  Dict("n_removed_params" => n_removed_params,
                "error_counter" => results_array[i][5],
                "param_subset" => string(param_subsets[i]),
                "dropped_fig" => results_array[i][1],
                "train_loss" => results_array[i][2],
                "test_loss" => results_array[i][3],
                "params" => results_array[i][4])

            results_df = push!(results_df, row_dict, promote=true)  
        end
       
        df_results = results_df

    #     max_length = maximum(length.(results_df.params))
    #     results_params_df = DataFrame([getindex.(results_df.params, i) for i in 1:max_length], :auto)
    #     rename!(results_params_df, param_names)

    #     df_results = hcat(select(results_df, Not(:params)), results_params_df)    

        # Filter subsets with NaN train or test loss values
        filter!(row -> !isnan(row.train_loss), df_results)
        filter!(row -> !isnan(row.test_loss), df_results)
     
        # Filter DataFrame df_result to only include rows that have train_loss lower than 1.1*minimum(train_loss)
        filter!(row -> row.train_loss < 1.1 * minimum(df_results.train_loss), df_results)

        # Sort DataFrame df_result by train_loss from low to high
        sort!(df_results, :train_loss)

        # Take the n_best_models param_subsets for this n_removed_params
        best_param_subsets_for_this_n_removed_params = eval.(Meta.parse.(df_results.param_subset))

        # # Add the best param_subset for this n_removed_params to the list of best_param_subsets for later test loss calculation
        # push!(best_param_subsets, best_param_subsets_for_this_n_removed_params...)

        # add the best params row for later use:     
        df_results_min = DataFrame(df_results[argmin(df_results.train_loss),:]) 
        df_results_all = vcat(df_results_all, df_results_min)
        
        # Update starting_param_subsets to be the best param_subsets for this n_removed_params
        global starting_param_subsets = values.(best_param_subsets_for_this_n_removed_params)

    end

    max_length = maximum(length.(df_results_all.params))
    results_params_df = DataFrame([getindex.(df_results_all.params, i) for i in 1:max_length], :auto)
    rename!(results_params_df, param_names)

    df_results_all = hcat(select(df_results_all, Not(:params)), results_params_df)   


    return df_results_all
end


function find_best_subset(data,
     vals_all_param_subsets,
     param_subsets_per_n_params,
     config,
     dir_path; 
     n_iter = 20,
     maxiter = 50_000,
     pct_bound = 0, 
     save_results = false, 
     print_prog = false, 
     all_parallel_subsets = false,
     )

    
    n_non_params = length(config["non_opt_params"])

     suffix = ""
    if all_parallel_subsets==true
        @info "starting all parallel subsets"
        # find the best subset parallel computation over all subsets (filtering out only na/inf loss)
        all_param_subsets_filter = filter_out_subsets_nan_values(data, vals_all_param_subsets, config) 
        @info "n param subsets after filtering: $(length(all_param_subsets_filter))"
        results_combined_df =  find_best_models_per_fig_all_subsets(data, all_param_subsets_filter, config, n_non_params=n_non_params, 
        n_iter = n_iter, maxiter = maxiter, save_results = save_results)
        suffix = "parallel"
        @info "finished finding best models per fig for all subsets, nrows = $(nrow(results_combined_df))"
    else
        @info "starting all subsets with filtering"
        # performing Denis's filtering 
        # find best model(s) for each figure
        figs = unique(data.source) 
        results_figs_dfs = pmap(
        dropped_fig -> find_best_subsets_per_fig(dropped_fig,param_subsets_per_n_params,data, config;
            n_iter =n_iter , maxiter = maxiter, print_prog = print_prog),
            figs,
        )
        results_combined_df = vcat(results_figs_dfs...)
        suffix = "filtering_subsets"
        @info "finished finding best models per fig"
    end

    if save_results==true
        save_csv_with_subdir(results_combined_df,
            joinpath(dir_path, "Cluster_results/$date_time"),    
            "PKM2_cv_best_model_per_figure_and_n_removed_param_$suffix.csv"
        )
    end

    ### find best number of parameters:x
    ## pct_bound - percent from minimal test loss 
    ## first: keep number of parameters with cv test loss < (1+ pct_bound/1000) * minimal cv test loss
    ## second: among them, choose the maximum n dropped params 
    best_n_removed_params = find_best_n_removed_params(results_combined_df,pct_bound= pct_bound, print_prog = print_prog)

    ### train all models with the best number using training data = all data and choose the best one with minimal train loss
    ## also save the data frame with all models
    best_param_subset_row = train_and_choose_best_subset(data, param_subsets_per_n_params, best_n_removed_params, config,dir_path, n_iter = n_iter, 
        maxiter = maxiter, print_prog = print_prog, save_results = save_results)

    #     max_length = maximum(length.(results_df.params))
    #     results_params_df = DataFrame([getindex.(results_df.params, i) for i in 1:max_length], :auto)
    #     rename!(results_params_df, param_names)
    #     results_df = hcat(select(results_df, Not(:params)), results_params_df)  


    if save_results == true
        save_csv_with_subdir(best_param_subset_row,
        joinpath(dir_path, "Cluster_results/$date_time"),     
         "PKM2_best_subset_$suffix.csv")
    end 

    return best_param_subset_row
end


function find_best_n_removed_params(df_results; pct_bound = 0, print_prog = false)
    # for each n_removed_params and dropped fig: keep the best model (minimal training_loss)
    function min_loss_row(group)
        return group[argmin(group.train_loss), :]
    end

    # Combine by 'n_removed_params' and 'dropped_fig' and apply the min_loss_row function
    df_results = combine(groupby(df_results, [:n_removed_params, :dropped_fig]), min_loss_row)

    # Calculate average test loss for each n_removed_parameters
    avg_values = combine(groupby(df_results, :n_removed_params), :test_loss => mean => :avg_test_loss)

    n_removed_params_almost_best = filter(row -> row.avg_test_loss <= (1+(pct_bound/100)) * minimum(avg_values.avg_test_loss), avg_values)

    if print_prog == true
        println("Avg CV error for each n removed params:")
        println(sort(avg_values, :avg_test_loss))
    end


    almost_best_with_maximum_removed_params = argmax(n_removed_params_almost_best.avg_test_loss)
    best_n_removed_params = n_removed_params_almost_best[almost_best_with_maximum_removed_params, :].n_removed_params

    if print_prog == true
        println("best n_removed_parameters: ", best_n_removed_params, ", PCT boundary: $pct_bound%")
    end

    return best_n_removed_params
end

function train_and_choose_best_subset(data,param_subsets_per_n_params,  best_n_removed_params, config, dir_path; n_iter = 1, maxiter = 50_000, save_results = false, print_prog = false)
    param_subsets = param_subsets_per_n_params[best_n_removed_params]
    #     param_subsets = [NamedTuple{Tuple(nt_param_choice_names)}(x) for x in unique(param_subsets)]
    param_subsets = filter_out_subsets_nan_values(data, param_subsets, config)

    results_array = pmap(
        param_subset -> train(data, config; n_iter = n_iter, nt_param_choice = param_subset, maxiter = maxiter),
        param_subsets,
    )


    results_df = DataFrame(
        n_removed_params = Int[], 
        error_counter = Int[],
        param_subset = String[], 
        train_loss = String[], 
        params = Vector{Float64}(undef,0)
    )


    for i in eachindex(results_array)
        row_dict =  Dict("n_removed_params" => best_n_removed_params,
            "error_counter" => results_array[i][3],
            "param_subset" => string(param_subsets[i]),
            "train_loss" => results_array[i][1],
            "params" => results_array[i][2] # before rescaling
    #            "params" => param_subset_select(param_rescaling(results_array[i][2])  , param_subsets[i])
        )

        results_df = push!(results_df, row_dict, promote=true)  
    end


 #     max_length = maximum(length.(results_df.params))
 #     results_params_df = DataFrame([getindex.(results_df.params, i) for i in 1:max_length], :auto)
 #     rename!(results_params_df, param_names)

 #     results_df = hcat(select(results_df, Not(:params)), results_params_df)   

    if save_results == true
        save_csv_with_subdir(results_df,
        joinpath(dir_path, "Cluster_results/$date_time"),   
        "PKM2_training_all_models_w_best_n_removed_params.csv")
    end

    best_param_subset = DataFrame(results_df[argmin(results_df.train_loss),:])
    if print_prog == true
        println("Best subset: $(best_param_subset.param_subset)")
    end

    return best_param_subset
end




function variable_selection(config, data, dir_path)
    # Format and print the time
    starting_time = now()
    formatted_time = Dates.format(starting_time, "HH:MM:SS")
    @info "starting time: $formatted_time"
     figs = unique(data.source) 
 
     # param subsets
     vals_all_param_subsets = create_param_subsets(config["param_constraints"], config["nt_param_choice_names"])
     # (0, 1, 2, 1, 1, 0, 0, 0, 0, 0) => (0, 0, 1, 1, 0)
     # number of terms at end of vals_all_param_subsets that are a3 where non-negative value != removing parameter
     n_non_params = length(config["non_opt_params"])
 
     # keep for each n removed params: all the subsets with this n:
     param_subsets_tuple = [(sum(values(param_subset)[1:end-n_non_params] .> 0), values(param_subset)) 
         for param_subset in vals_all_param_subsets]
 
     param_subsets_per_n_params = Dict{Int, Vector}()
 
     for (key, value) in param_subsets_tuple
         if haskey(param_subsets_per_n_params, key)
             push!(param_subsets_per_n_params[key], value)
         else
             param_subsets_per_n_params[key] = [value]
         end
     end
     # Dict{Int64, Vector} with 8 entries:
     #   0 => [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 1, 0), (0, 0, …
     #   4 => [(0, 1, 1, 1, 1, 0, 0, 0, 0, 0), (0, 1, 2, 1, 1, 0, 0, 0, 0, 0), (0, 1, …
     #   5 => [(0, 1, 1, 1, 1, 1, 0, 0, 0, 0), (0, 1, 2, 1, 1, 1, 0, 0, 0, 0), (0, 1, …
     #   6 => [(0, 1, 1, 1, 1, 1, 1, 0, 0, 0), (0, 1, 2, 1, 1, 1, 1, 0, 0, 0), (0, 1, …
     #   2 => [(0, 1, 1, 0, 0, 0, 0, 0, 0, 0), (0, 1, 2, 0, 0, 0, 0, 0, 0, 0), (0, 1, …
     #   7 => [(0, 1, 1, 1, 1, 1, 1, 1, 0, 0), (0, 1, 2, 1, 1, 1, 1, 1, 0, 0), (0, 1, …
     #   3 => [(0, 1, 1, 1, 0, 0, 0, 0, 0, 0), (0, 1, 2, 1, 0, 0, 0, 0, 0, 0), (0, 1, …
     #   1 => [(0, 1, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0, 0, 0, 0, 0), (0, 0, …

    #  best_param_subset_row = @time find_best_subset(data,vals_all_param_subsets, param_subsets_per_n_params, config,dir_path,  n_iter = config["fitting_params"]["n_iter"], maxiter = config["fitting_params"]["maxiter"], 
    #      print_prog = true, save_results = true,  all_parallel_subsets = false)

    best_param_subset_row = @time find_best_subset(data,vals_all_param_subsets, param_subsets_per_n_params, config,dir_path,  n_iter = config["fitting_params"]["n_iter"], maxiter = config["fitting_params"]["maxiter"], 
         print_prog = true, save_results = true,  all_parallel_subsets = true)

     # Format and print the time
     ending_time = now()
     formatted_time = Dates.format(ending_time, "HH:MM:SS")
     @info "ending time: $formatted_time"
 
     duration_minutes = (ending_time - starting_time) / Dates.Minute(1)
     @info "during time: $duration_minutes"
end
