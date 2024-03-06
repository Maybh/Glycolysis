using CMAEvolutionStrategy, JSON, Statistics, Distributions, LinearAlgebra, MixedModels

include("Utils.jl")
include("Enzyme.jl")

function loss_likelihood(
    kinetic_params,
    rate_data,
    fig_point_indexes, 
    config;
    nt_param_choice = nothing,
    optimization_param_names = nothing
    )

    nt_param_choice_names = config["nt_param_choice_names"]
    if isnothing(nt_param_choice)
        nt_param_choice = (; zip(nt_param_choice_names, zeros(Int64, length(nt_param_choice_names)))...)
    end
    
    if isnothing(optimization_param_names)
        kinetic_params = NamedTuple{Tuple(config["param_names"])}(kinetic_params)
    # if optimization_param_names isn't nothing, it means we have partial list of the parameters 
    else
        kinetic_params = NamedTuple{Tuple(Symbol.(optimization_param_names))}(kinetic_params)
    end 
    
    # rescaling and param_constraints
    kinetic_params = param_rescaling_from_conf(kinetic_params, config["rescaling"])
    kinetic_params = param_subset_select_from_conf(kinetic_params, nt_param_choice, 
        config["param_constraints"], config["param_names"], config["constant_params"])
    
    log_actual_vs_pred = Vector{Float64}(undef, length(rate_data.Rate))
    for i = 1:length(rate_data.Rate)
        log_actual_vs_pred[i] = log(
                rate_data.Rate[i] / rate_PKM2(
                    rate_data.PEP[i],
                    rate_data.ADP[i],
                    rate_data.Pyruvate[i],
                    rate_data.ATP[i],
                    rate_data.F16BP[i],
                    rate_data.Phenylalanine[i],
                    kinetic_params,
                ) ,
        )
    end
    
    df = DataFrame(figure=rate_data.fig_num, log_ratios=log_actual_vs_pred)
    formula = @formula(log_ratios ~ (1|figure))
    model = fit!(LinearMixedModel(formula, df))
    
    return -loglikelihood(model) # Negative log-likelihood for minimization
end

"Loss function used for fitting"
function loss_likelihood_old(
    kinetic_params,
    rate_data,
    fig_point_indexes, 
    config;
    nt_param_choice = nothing,
    optimization_param_names = nothing
    )

    nt_param_choice_names = config["nt_param_choice_names"]
    if isnothing(nt_param_choice)
        nt_param_choice = (; zip(nt_param_choice_names, zeros(Int64, length(nt_param_choice_names)))...)
    end
    
    if isnothing(optimization_param_names)
        kinetic_params = NamedTuple{Tuple(config["param_names"])}(kinetic_params)
    # if optimization_param_names isn't nothing, it means we have partial list of the parameters 
    else
        kinetic_params = NamedTuple{Tuple(Symbol.(optimization_param_names))}(kinetic_params)
    end 
    
    # rescaling and param_constraints
    kinetic_params = param_rescaling_from_conf(kinetic_params, config["rescaling"])
    kinetic_params = param_subset_select_from_conf(kinetic_params, nt_param_choice, 
        config["param_constraints"], config["param_names"], config["constant_params"])
    
    log_likelihood_total = 0.0
    for j in 1:maximum(rate_data.fig_num)
        fig_indexes = fig_point_indexes[j]
        figure_data = NamedTuple{keys(rate_data)}(tuple([getfield(rate_data, k)[fig_indexes] for k in keys(rate_data)]...))
        n_j = length(figure_data.Rate)
        

        # Construct the covariance matrix
        Sigma_j = kinetic_params.sigma_alpha * ones(n_j, n_j) + kinetic_params.sigma * I(n_j)
       
        log_ratios = Vector{Float64}(undef, length(figure_data.Rate))
        for i = 1:n_j
            log_ratios[i] = log(
                figure_data.Rate[i] / rate_PKM2(
                    figure_data.PEP[i],
                    figure_data.ADP[i],
                    figure_data.Pyruvate[i],
                    figure_data.ATP[i],
                    figure_data.F16BP[i],
                    figure_data.Phenylalanine[i],
                    kinetic_params,
                ) ,
            )
        end
        
        # Calculate the multivariate normal log-likelihood for figure j
        likelihood_figure = MvNormal(zeros(n_j), Sigma_j)
        log_likelihood_total += logpdf(likelihood_figure, log_ratios)
    end
    return -log_likelihood_total  # Negative log-likelihood for minimization
end

function loss_sse(
        kinetic_params,
        rate_data,
        fig_point_indexes, 
        config;
        nt_param_choice = nothing,
        optimization_param_names = nothing, 
            
    )

    # @info "kinetic params: $kinetic_params"
    # @info "nt_param_choice: $nt_param_choice"
    # @info "optimization_param_names: $optimization_param_names"

    nt_param_choice_names = config["nt_param_choice_names"]
    if isnothing(nt_param_choice)
        nt_param_choice = (; zip(nt_param_choice_names, zeros(Int64, length(nt_param_choice_names)))...)
    end
  
    ### Convert kinetic params to NamedTuple:
    # if optimization_param_names is nothing use param_names and create NamedTuple
    if isnothing(optimization_param_names)
        kinetic_params = NamedTuple{Tuple(config["param_names"])}(kinetic_params)
    # if optimization_param_names isn't nothing, it means we have partial list of the parameters 
    else
        kinetic_params = NamedTuple{Tuple(Symbol.(optimization_param_names))}(kinetic_params)
    end 
 
    # rescaling and param_constraints
    kinetic_params = param_rescaling_from_conf(kinetic_params, config["rescaling"])
    kinetic_params = param_subset_select_from_conf(kinetic_params, nt_param_choice, 
            config["param_names"], config["constant_params"])

    loss = zero(eltype(kinetic_params))
    log_pred_vs_data_ratios = Vector{Float64}(undef, length(rate_data.Rate))
    #precalculate rate_PKM2() for all points as it is expensive and reuse it for weights and loss
    for i = 1:length(rate_data.Rate)
        rate_data_i = extract_data_at_index(rate_data, i)
        log_pred_vs_data_ratios[i] = log(
            rate_function(rate_data_i, 
            config["substrates"], 
            config["products"], 
            config["regulators"],
            config["reg_binding_sites"],
            config["n"],
            kinetic_params,
            ) / rate_data.Rate[i],
        )
    end
    #Calculate figures weights and loss on per figure basis
    for i = 1:maximum(rate_data.fig_num)
        # calculate Vmax weights for each figure which have analytical solution as ratio of gemetric means of data vs prediction
        log_weight = mean(-log_pred_vs_data_ratios[fig_point_indexes[i]])
        loss += sum(abs2.(log_weight .+ log_pred_vs_data_ratios[fig_point_indexes[i]]))
    end

    return loss / length(rate_data.Rate)
end


function loss( 
    kinetic_params,
    rate_data,
    fig_point_indexes, 
    config;
    nt_param_choice = nothing,
    optimization_param_names = nothing
    )
 
    if config["loss_choice"] == "sse"
        loss_res = loss_sse(
            kinetic_params,
            rate_data,
            fig_point_indexes, 
            config;
            nt_param_choice = nt_param_choice,
            optimization_param_names = optimization_param_names
        )
    elseif config["loss_choice"] == "likelihood"
        loss_res = loss_likelihood(
            kinetic_params,
            rate_data,
            fig_point_indexes, 
            config;
            nt_param_choice = nt_param_choice,
            optimization_param_names = optimization_param_names
        )
    end
    return loss_res
end

function get_loss_function(loss_choice = "sse")
    if loss_choice == "sse"
        return loss_sse
    elseif loss_choice == "likelihood"
        return loss_likelihood
    end
end

"function to calculate train loss without a figure and test loss on removed figure"
function fig_loocv(
        fig,
        data,
        config;
        n_iter = 20,
        maxiter = 50_000,
        nt_param_choice = nothing,
    )
 
    nt_param_choice_names = config["nt_param_choice_names"]
    if isnothing(nt_param_choice)
        nt_param_choice = (; zip(nt_param_choice_names, zeros(Int64, length(nt_param_choice_names)))...)
    end
   
    # Drop selected figure from data
    train_data = data[data.source.!=fig, :]
    test_data = data[data.source.==fig, :]
    # Calculate fit
    train_res = train(train_data, config; n_iter = n_iter, nt_param_choice = nt_param_choice, maxiter = maxiter)
    test_res = test(test_data, train_res[2], config; nt_param_choice = nt_param_choice, 
                        optimization_param_names = train_res[4])

    ## return params after rescaling and param subset select: 
    est_params = NamedTuple{Tuple(train_res[4])}(train_res[2])
    est_params = param_rescaling_from_conf(est_params, config["rescaling"])
    est_params = param_subset_select_from_conf(est_params, nt_param_choice, config["param_constraints"],
        config["param_names"], config["constant_params"])
    # for now keep only list of the values. later I want to fix it to send Named Tuple (and change code for that)
    est_params_lst = collect(values(est_params))

    return fig, train_res[1], test_res, est_params_lst, train_res[3]
end

function train(
        data,
        config;
        n_iter = 20,
        nt_param_choice = nothing,
        maxiter = 50_000
    )

    param_names = config["param_names"]
    metab_names =  Symbol.(collect(keys(config["metab_names"])))
    nt_param_choice_names = config["nt_param_choice_names"]
    
    if isnothing(nt_param_choice)
        nt_param_choice = (; zip(nt_param_choice_names, zeros(Int64, length(nt_param_choice_names)))...)
    end
   
    # Add a new column to data to assign an integer to each source/figure from publication
    data.fig_num = vcat(
        [
            i * ones(Int64, count(==(unique(data.source)[i]), data.source)) for
            i = 1:length(unique(data.source))
        ]...,
    )
    # Convert DF to NamedTuple for better type stability / speed
    rate_data = Tables.columntable(data[.!isnan.(data.Rate), [:Rate, metab_names..., :fig_num]])
    # Make a vector containing indexes of points corresponding to each figure
    fig_point_indexes = [findall(data.fig_num .== i) for i in unique(data.fig_num)]

    loss = get_loss_function(config["loss_choice"])
    
    # Check if nt_param_choice makes loss returns NaN and abort early if it does. The latter could happens due to nt_param_choice making params=Inf
    if isnan(
        loss(
            5 .* ones(length(param_names)),
            rate_data,
            fig_point_indexes, 
            config,
            nt_param_choice = nt_param_choice,
        ),
    )
        println("Loss returns NaN for this param combo")
        return Inf, fill(NaN, length(param_names)), n_iter, nothing
    end
  
    keys_to_exclude = config["non_opt_params"]
    # n_params_for_optimization = length(param_names) - sum((nt_param_choice[key] > 0) for key in fieldnames(typeof(nt_param_choice)) if !(key in keys_to_exclude))
    param_names_for_optimization = get_param_names_for_optimization(nt_param_choice, config["param_constraints"], 
        config["param_names"], config["constant_params"], config["non_opt_params"] )
    n_params_for_optimization = length(param_names_for_optimization)


    # s = Sobol.SobolSeq([0.0 for i = 1:length(param_names)], [10.0 for i = 1:length(param_names)])
    # Sobol.skip(s, n_iter)
    solns = []
    error_counter = 0
    for i = 1:n_iter
        # x0 = Sobol.next!(s)
        x0 = 10 .* rand(n_params_for_optimization)
        sol = try
            minimize(
                x -> loss(x, rate_data, fig_point_indexes, config; 
                                nt_param_choice = nt_param_choice, 
                                optimization_param_names = param_names_for_optimization),
                x0,
                0.01,
                lower = repeat([0.0], length(x0)),
                upper = repeat([10.0], length(x0)),
                popsize = 4 * (4 + floor(Int, 3 * log(length(x0)))),
                maxiter = maxiter,
                verbosity = 0,
                ftol = 1e-10,
            )
        catch errorbest_param_subset_row
            # bypass rare errors (~1 in 10,000 runs) where the minimize() fails to converge with "ArgumentError: matrix contains Infs or NaNs"
            if isa(error, ArgumentError)
                println(error)
                error_counter += 1
                sol = Inf
            end
            println(error)
        end
        push!(solns, sol)
    end

    solns = [sol for sol in solns if !isnothing(sol) && !isinf(sol) && !isnan(sol)]
    if isempty(solns)
        println("All of the iterations of fits for this param combo return NaN or Inf")
        return Inf, fill(NaN, n_params_for_optimization), error_counter, param_names_for_optimization
    end
    index_best_sol = argmin([fbest(sol) for sol in solns])

    best_sol = try
        minimize(
            x -> loss(x, rate_data, fig_point_indexes, config; 
                            nt_param_choice = nt_param_choice, 
                            optimization_param_names = param_names_for_optimization),
            xbest(solns[index_best_sol]),
            0.001,
            lower = repeat([0.0], length(xbest(solns[index_best_sol]))),
            upper = repeat([10.0], length(xbest(solns[index_best_sol]))),
            popsize = 4 * (4 + floor(Int, 3 * log(length(xbest(solns[index_best_sol]))))),
            maxiter = maxiter,
            verbosity = 0,
            ftol = 1e-14,
        )
    catch error
        # bypass rare errors where the minimize() fails to converge with "ArgumentError: matrix contains Infs or NaNs"
        if isa(error, ArgumentError)
            println(error)
            error_counter += 1
            best_sol = solns[index_best_sol]
        end
    end
    return fbest(best_sol), xbest(best_sol), error_counter, param_names_for_optimization
end

    

"function to test PKM2 model fit to data"
function test(
        data,
        kinetic_params, 
        config;
        nt_param_choice = nothing,
        optimization_param_names = nothing
    )
        
    nt_param_choice_names = config["nt_param_choice_names"]
    
    if isnothing(nt_param_choice)
        nt_param_choice = (; zip(nt_param_choice_names, zeros(Int64, length(nt_param_choice_names)))...)
    end
    param_names = config["param_names"]
    metab_names = config["metab_names"]

    # Add a new column to data to assign an integer to each source/figure from publication
    data.fig_num = vcat(
        [
            i * ones(Int64, count(==(unique(data.source)[i]), data.source)) for
            i = 1:length(unique(data.source))
        ]...,
    )

    # Convert DF to NamedTuple for better type stability / speed
    rate_data = Tables.columntable(data[.!isnan.(data.Rate), [:Rate, metab_names..., :fig_num]])

    # Make a vector containing indexes of points corresponding to each figure
    fig_point_indexes = [findall(data.fig_num .== i) for i in unique(data.fig_num)]

    # Check if one of the parameters is NaN
    if any(isnan.(kinetic_params))
        println("One of the kinetic parameters is NaN")
        return Inf
    end
    
    test_loss = loss(kinetic_params, rate_data, fig_point_indexes, config;
    nt_param_choice = nt_param_choice,
     optimization_param_names = optimization_param_names)
    
    if isnan(test_loss)
        println("Loss returns NaN for this param combo")
        return Inf
    end

    return test_loss
end

"get all the parameter names for optimization"
function get_param_names_for_optimization(nt_param_choice, param_constraints, param_names, constant_params, non_opt_params)
    # get keys name
    constant_param_names = collect(keys(constant_params))
    # get constrains parameters name
    constraints_param_names = vcat([collect(keys(param_constraints[string(name)][string(value)])) for (name, value) in pairs(nt_param_choice) if (haskey(param_constraints, string(name)) & value > 0)]...)

    pramas_to_ignore = vcat(constant_param_names, constraints_param_names, non_opt_params)
    
    # Remove the intersection of pramas_to_ignore from param_names
    param_names_for_optimization = setdiff(param_names, Symbol.(pramas_to_ignore))
    return param_names_for_optimization
end

# Implement isnan and isinf for CMAEvolutionStrategy.Optimizer type
Base.isnan(sol::CMAEvolutionStrategy.Optimizer) = isnan(fbest(sol))
Base.isinf(sol::CMAEvolutionStrategy.Optimizer) = isinf(fbest(sol))