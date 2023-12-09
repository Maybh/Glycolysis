using CMAEvolutionStrategy

#=
This files contains most of the functions necessary to identify the best MWC rate 
equation for PKM2 and fit it to data.

To adopt this to another MWC enzyme change the following functions:
    rate_PKM2() and all of the places where it is called
    param_names
    metab_names
    param_subset_select()
    nt_param_choice_names
    and maybe param_rescaling()
=#


"Names of parameters. Make sure it matches exactly allocation of p in rate_PKM2()"
param_names = [
    :L,
    :Vmax_a,
    :Vmax_i,
    :K_a_PEP,
    :K_a_ADP,
    :K_a_Pyruvate,
    :K_a_ATP,
    :K_a_F16BP,
    :K_a_Phenylalanine,
    :K_i_PEP,
    :K_i_ADP,
    :K_i_Pyruvate,
    :K_i_ATP,
    :K_i_F16BP,
    :K_i_Phenylalanine,
    :alpha_PEP_ATP,
    :alpha_Pyruvate_ADP,
]

"Names of PKM substrate, products and regulators"
metab_names = [:PEP, :ADP, :Pyruvate, :ATP, :F16BP, :Phenylalanine]

"function to calculate rate of PKM2"
function rate_PKM2(PEP, ADP, Pyruvate, ATP, F16BP, Phenylalanine, p)
    L,
    Vmax_a,
    Vmax_i,
    K_a_PEP,
    K_a_ADP,
    K_a_Pyruvate,
    K_a_ATP,
    K_a_F16BP,
    K_a_Phenylalanine,
    K_i_PEP,
    K_i_ADP,
    K_i_Pyruvate,
    K_i_ATP,
    K_i_F16BP,
    K_i_Phenylalanine,
    alpha_PEP_ATP,
    alpha_Pyruvate_ADP = p

    Keq = 20_000
    Vmax_a = 1.0


    # Vmax_i = Vmax_a
    # K_i_Pyruvate = K_a_Pyruvate
    # K_i_ATP = K_a_ATP
    # K_i_ADP = K_a_ADP
    # K_a_Phenylalanine = Inf
    # K_i_F16BP = Inf
    # alpha_PEP_ATP = 0.0
    # alpha_Pyruvate_ADP = 1.0

    Vmax_a_rev =
        (K_a_ATP != Inf && K_a_Pyruvate != Inf) ?
        Vmax_a * K_a_ATP * K_a_Pyruvate / (Keq * K_a_PEP * K_a_ADP) : 0.0
    Vmax_i_rev =
        (K_i_ATP != Inf && K_i_Pyruvate != Inf) ?
        Vmax_i * K_i_ATP * K_i_Pyruvate / (Keq * K_i_PEP * K_i_ADP) : 0.0

    Z_a_cat = (
        1 +
        (PEP / K_a_PEP) +
        (ATP / K_a_ATP) +
        (Pyruvate / K_a_Pyruvate) +
        (ADP / K_a_ADP) +
        (PEP / K_a_PEP) * (ADP / K_a_ADP) +
        (Pyruvate / K_a_Pyruvate) * (ATP / K_a_ATP) +
        alpha_PEP_ATP * (PEP / K_a_PEP) * (ATP / K_a_ATP) +
        alpha_Pyruvate_ADP * (Pyruvate / K_a_Pyruvate) * (ADP / K_a_ADP)
    )
    Z_i_cat = (
        1 +
        (PEP / K_i_PEP) +
        (ATP / K_i_ATP) +
        (Pyruvate / K_i_Pyruvate) +
        (ADP / K_i_ADP) +
        (PEP / K_i_PEP) * (ADP / K_i_ADP) +
        (Pyruvate / K_i_Pyruvate) * (ATP / K_i_ATP) +
        alpha_PEP_ATP * (PEP / K_i_PEP) * (ATP / K_i_ATP) +
        alpha_Pyruvate_ADP * (Pyruvate / K_i_Pyruvate) * (ADP / K_i_ADP)
    )
    Z_a_reg = ((1 + F16BP / K_a_F16BP) * (1 + Phenylalanine / K_a_Phenylalanine))
    Z_i_reg = ((1 + F16BP / K_i_F16BP) * (1 + Phenylalanine / K_i_Phenylalanine))

    Rate =
        (
            (
                Vmax_a * (PEP / K_a_PEP) * (ADP / K_a_ADP) -
                Vmax_a_rev * (Pyruvate / K_a_Pyruvate) * (ATP / K_a_ATP)
            ) *
            (Z_a_cat^3) *
            (Z_a_reg^4) +
            L *
            (
                Vmax_i * (PEP / K_i_PEP) * (ADP / K_i_ADP) -
                Vmax_i_rev * (Pyruvate / K_i_Pyruvate) * (ATP / K_i_ATP)
            ) *
            (Z_i_cat^3) *
            (Z_i_reg^4)
        ) / ((Z_a_cat^4) * (Z_a_reg^4) + L * (Z_i_cat^4) * (Z_i_reg^4))

    return Rate
end

"function to calculate train loss without a figure and test loss on removed figure"
function fig_loocv_PKM2(
    fig,
    data,
    param_names;
    n_iter = 1,
    nt_param_choice = (; zip(nt_param_choice_names, zeros(Int64, length(nt_param_choice_names)))...),
)
    # Drop selected figure from data
    train_data = data[data.source.!=fig, :]
    test_data = data[data.source.==fig, :]
    # Calculate fit
    train_res = train_PKM2(train_data, param_names; n_iter = n_iter, nt_param_choice = nt_param_choice)
    test_res = test_PKM2(test_data, train_res[2]; nt_param_choice = nt_param_choice)
    return fig, train_res[1], test_res, train_res[2], train_res[3]
end

"function to fit PKM2 model to data"
function train_PKM2(
    data,
    param_names;
    n_iter = 20,
    nt_param_choice = (; zip(nt_param_choice_names, zeros(Int64, length(nt_param_choice_names)))...),
)
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

    # Check if nt_param_choice makes loss returns NaN and abort early if it does. The latter could happens due to nt_param_choice making params=Inf
    if isnan(
        loss_PKM2!(
            5 .* ones(length(param_names)),
            rate_data,
            fig_point_indexes;
            nt_param_choice = nt_param_choice,
        ),
    )
        println("Loss returns NaN for this param combo")
        return Inf, fill(NaN, length(param_names)), n_iter
    end

    # s = Sobol.SobolSeq([0.0 for i = 1:length(param_names)], [10.0 for i = 1:length(param_names)])
    # Sobol.skip(s, n_iter)
    solns = []
    error_counter = 0
    for i = 1:n_iter
        # x0 = Sobol.next!(s)
        x0 = 10 .* rand(length(param_names))
        sol = try
            minimize(
                x -> loss_PKM2!(x, rate_data, fig_point_indexes; nt_param_choice = nt_param_choice),
                x0,
                0.01,
                lower = repeat([0.0], length(x0)),
                upper = repeat([10.0], length(x0)),
                popsize = 4 * (4 + floor(Int, 3 * log(length(x0)))),
                maxiter = 100,
                verbosity = 0,
                ftol = 1e-10,
            )
        catch error
            # bypass rare errors (~1 in 10,000 runs) where the minimize() fails to converge with "ArgumentError: matrix contains Infs or NaNs"
            if isa(error, ArgumentError)
                println(error)
                error_counter += 1
                sol = Inf
            end
        end
        push!(solns, sol)
    end

    solns = [sol for sol in solns if !isinf(sol) && !isnan(sol)]
    if isempty(solns)
        println("All of the iterations of fits for this param combo return NaN or Inf")
        return Inf, fill(NaN, length(param_names)), error_counter
    end
    index_best_sol = argmin([fbest(sol) for sol in solns])
    best_sol = try
        minimize(
            x -> loss_PKM2!(x, rate_data, fig_point_indexes; nt_param_choice = nt_param_choice),
            xbest(solns[index_best_sol]),
            0.001,
            lower = repeat([0.0], length(xbest(solns[index_best_sol]))),
            upper = repeat([10.0], length(xbest(solns[index_best_sol]))),
            popsize = 4 * (4 + floor(Int, 3 * log(length(xbest(solns[index_best_sol]))))),
            maxiter = 100,
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
    return fbest(best_sol), xbest(best_sol), error_counter
end

"function to test PKM2 model fit to data"
function test_PKM2(
    data,
    params;
    nt_param_choice = (; zip(nt_param_choice_names, zeros(Int64, length(nt_param_choice_names)))...),
)
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

    #deepcopy params so they are not modified in place
    kinetic_params = deepcopy(params)

    # Check if nt_param_choice makes loss returns NaN and abort early if it does. The latter could happens due to nt_param_choice making params=Inf
    if isnan(
        loss_PKM2!(
            5 .* ones(length(param_names)),
            rate_data,
            fig_point_indexes;
            nt_param_choice = nt_param_choice,
        ),
    )
        println("Loss returns NaN for this param combo")
        return Inf
    end

    # Check if one of the parameters is NaN
    if any(isnan.(kinetic_params))
        println("One of the kinetic parameters is NaN")
        return Inf
    end

    test_loss = loss_PKM2!(kinetic_params, rate_data, fig_point_indexes; nt_param_choice = nt_param_choice)
    return test_loss
end

"Loss function used for fitting"
function loss_PKM2!(
    kinetic_params,
    rate_data,
    fig_point_indexes;
    nt_param_choice = (; zip(nt_param_choice_names, zeros(Int64, length(nt_param_choice_names)))...),
)
    kinetic_params .= param_rescaling(kinetic_params)
    kinetic_params .= param_subset_select(kinetic_params, nt_param_choice)
    loss = zero(eltype(kinetic_params))
    log_pred_vs_data_ratios = Vector{Float64}(undef, length(rate_data.Rate))

    #precalculate rate_PKM2() for all points as it is expensive and reuse it for weights and loss
    for i = 1:length(rate_data.Rate)
        log_pred_vs_data_ratios[i] = log(
            rate_PKM2(
                rate_data.PEP[i],
                rate_data.ADP[i],
                rate_data.Pyruvate[i],
                rate_data.ATP[i],
                rate_data.F16BP[i],
                rate_data.Phenylalanine[i],
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

"Function to calculate weight for each figure based on geometric mean of data vs prediction for plotting"
function figure_weights_PKM2(kinetic_params, rate_data, fig_point_indexes)
    nt_rate_data = Tables.columntable(rate_data[:, [:Rate, metab_names..., :fig_num]])

    #precalculate rate_PKM2() for all points as it is expensive and reuse it for weights
    log_pred_vs_data_ratios = Vector{Float64}(undef, length(nt_rate_data.Rate))
    for i = 1:length(nt_rate_data.Rate)
        log_pred_vs_data_ratios[i] = log(
            rate_PKM2(
                nt_rate_data.PEP[i],
                nt_rate_data.ADP[i],
                nt_rate_data.Pyruvate[i],
                nt_rate_data.ATP[i],
                nt_rate_data.F16BP[i],
                nt_rate_data.Phenylalanine[i],
                kinetic_params,
            ) / nt_rate_data.Rate[i],
        )
    end
    log_weights = Vector{Float64}(undef, maximum(nt_rate_data.fig_num))
    for i = 1:maximum(nt_rate_data.fig_num)
        # calculate Vmax weights for each figure which have analytical solution as ratio of gemetric means of data vs prediction
        log_weights[i] = mean(-log_pred_vs_data_ratios[fig_point_indexes[i]])
    end

    # Calculate a vector of weights for each point based on figure weights
    weights = Vector{Float64}(undef, length(nt_rate_data.Rate))
    for i = 1:length(nt_rate_data.Rate)
        weights[i] = exp.(log_weights[nt_rate_data.fig_num[i]])
    end

    return weights
end

"Rescaling function for parameter so that all fitting param fall on [0, 10] scale"
function param_rescaling(p)
    new_p = similar(p)
    for i in eachindex(p)
        if i == 1 # L
            new_p[i] = 10^(-5) * 10^(10 * p[i] / 10)
        elseif i >= 2 && i < 4 # V_a & V_i
            new_p[i] = 10^(-3) * 10^(3 * p[i] / 10)
        elseif i >= 4 && i < 16
            new_p[i] = 10^(-10) * 10^(13 * p[i] / 10)
        elseif i >= 16 #alphas
            p[i] >= 5.0 ? new_p[i] = 1.0 : new_p[i] = 0.0
        end
    end
    return new_p
end

"""
Function to convert parameter vector to vector where some params are equal to 0, Inf or each other based on nt_param_choice.
Make sure the name of params for p is exactly the same as for rate_PKM2() and param_names.
"""
function param_subset_select(p, nt_param_choice)
    L,
    Vmax_a,
    Vmax_i,
    K_a_PEP,
    K_a_ADP,
    K_a_Pyruvate,
    K_a_ATP,
    K_a_F16BP,
    K_a_Phenylalanine,
    K_i_PEP,
    K_i_ADP,
    K_i_Pyruvate,
    K_i_ATP,
    K_i_F16BP,
    K_i_Phenylalanine,
    alpha_PEP_ATP,
    alpha_Pyruvate_ADP = p

    Vmax_a = 1.0

    if nt_param_choice.L == 1
        L = 0.0
    end

    if nt_param_choice.Vmax == 1
        Vmax_i = Vmax_a
    end

    if nt_param_choice.K_PEP == 1
        K_i_PEP = K_a_PEP
    elseif nt_param_choice.K_PEP == 2
        K_a_PEP = Inf
    elseif nt_param_choice.K_PEP == 3
        K_i_PEP = Inf
    end

    if nt_param_choice.K_ADP == 1
        K_i_ADP = K_a_ADP
    elseif nt_param_choice.K_ADP == 2
        K_a_ADP = Inf
    elseif nt_param_choice.K_ADP == 3
        K_i_ADP = Inf
    end

    if nt_param_choice.K_Pyruvate == 1
        K_i_Pyruvate = K_a_Pyruvate
    elseif nt_param_choice.K_Pyruvate == 2
        K_a_Pyruvate = Inf
    elseif nt_param_choice.K_Pyruvate == 3
        K_i_Pyruvate = Inf
    end

    if nt_param_choice.K_ATP == 1
        K_i_ATP = K_a_ATP
    elseif nt_param_choice.K_ATP == 2
        K_a_ATP = Inf
    elseif nt_param_choice.K_ATP == 3
        K_i_ATP = Inf
    end

    if nt_param_choice.K_F16BP == 1
        K_i_F16BP = K_a_F16BP
    elseif nt_param_choice.K_F16BP == 2
        K_a_F16BP = Inf
    elseif nt_param_choice.K_F16BP == 3
        K_i_F16BP = Inf
    end

    if nt_param_choice.K_Phenylalanine == 1
        K_i_Phenylalanine = K_a_Phenylalanine
    elseif nt_param_choice.K_Phenylalanine == 2
        K_a_Phenylalanine = Inf
    elseif nt_param_choice.K_Phenylalanine == 3
        K_i_Phenylalanine = Inf
    end

    if nt_param_choice.alpha_PEP_ATP == 0
        alpha_PEP_ATP = 0.0
    elseif nt_param_choice.alpha_PEP_ATP == 1
        alpha_PEP_ATP = 1.0
    end

    if nt_param_choice.alpha_Pyruvate_ADP == 0
        alpha_Pyruvate_ADP = 0.0
    elseif nt_param_choice.alpha_Pyruvate_ADP == 1
        alpha_Pyruvate_ADP = 1.0
    end

    new_p = L,
    Vmax_a,
    Vmax_i,
    K_a_PEP,
    K_a_ADP,
    K_a_Pyruvate,
    K_a_ATP,
    K_a_F16BP,
    K_a_Phenylalanine,
    K_i_PEP,
    K_i_ADP,
    K_i_Pyruvate,
    K_i_ATP,
    K_i_F16BP,
    K_i_Phenylalanine,
    alpha_PEP_ATP,
    alpha_Pyruvate_ADP
    return new_p
end

nt_param_choice_names = [
    :L,
    :Vmax,
    :K_PEP,
    :K_ADP,
    :K_Pyruvate,
    :K_ATP,
    :K_F16BP,
    :K_Phenylalanine,
    :alpha_PEP_ATP,
    :alpha_Pyruvate_ADP,
]

# Implement isnan and isinf for CMAEvolutionStrategy.Optimizer type
Base.isnan(sol::CMAEvolutionStrategy.Optimizer) = isnan(fbest(sol))
Base.isinf(sol::CMAEvolutionStrategy.Optimizer) = isinf(fbest(sol))