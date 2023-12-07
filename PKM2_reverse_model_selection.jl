#cd(@__DIR__)
using Pkg
#Pkg.activate(".")
Pkg.instantiate()

using Distributed, ClusterManagers, Dates

#Increase time before worker terminates without hearing from master process. Useful for large number of cores.
ENV["JULIA_WORKER_TIMEOUT"] = 600.0

# ENV["SLURM_JOB_CPUS_PER_NODE"] give number of cores in "40(x4),32,20" format
# "40(x2),32,20" can be parsed to 132 using code below to get total number of allocated cores

subs = Dict("x" => "*", "(" => "", ")" => "");
np = sum(eval(Meta.parse(replace(ENV["SLURM_JOB_CPUS_PER_NODE"], r"x|\(|\)" => s -> subs[s]))))
addprocs(SlurmManager(np); exeflags = "--project")

#Uncomment code below to run on local machine or uncomment code above to run on cluster
#If run locally adjust number of workers to number of cores on your machine
# using Distributed, Dates
# addprocs(8; exeflags = "--project")

#= 
This code perform forward model selection to fit different version of rate_PKM2() equation to full data.
This should be run on cluster
=#

@everywhere using CMAEvolutionStrategy, Optim, Random, Sobol, DataFrames, Statistics, CSV

#Load function of PKM2 fitting
@everywhere include("PKM2_fitting_functions.jl")

#Load enzyme kinetic data data
PKM2_data_for_fit = CSV.read("PKM2_data/PKM2_data.csv", DataFrame)

#Add source column that uniquely identifies a figure from publication
PKM2_data_for_fit.source .= PKM2_data_for_fit.Article .* "_" .* PKM2_data_for_fit.Fig

# Make an array of all possible param combination
a1 = [0, 1]
a2 = [0, 1, 2, 3]
a3 = [0, 1]

# param subsets
vals_all_param_subsets = collect(Iterators.product(0, a1, a2, a2, a2, a2, a2, a2, a3, a3))
# number of terms at end of vals_all_param_subsets that are a3 where non-negative value != removing parameter
n_non_params = 2
starting_param_subsets = [x for x in vals_all_param_subsets if sum(values(x[1:end-n_non_params]) .> 0) == 7]

# Initialize array to save the best models from each round
best_param_subsets = NamedTuple[]

# Loop over number of parameters
for n_removed_params = 7:-1:0
    # Choose param_subsets with n_removed_params removed that also contain zero elements from starting_param_subsets
    #make a mask of zero elements from starting_param_subsets
    starting_param_subset_masks =
        unique([starting_param_subset .> 0 for starting_param_subset in starting_param_subsets])
    #select all subsets with n_removed_params non-zero elements
    subset_of_all_param_subsets_w_n_removed_params = [
        param_subsets for
        param_subsets in vals_all_param_subsets if sum((param_subsets .> 0)[1:end-n_non_params]) == n_removed_params
    ]
    #actually choose param_subsets with n_removed_params of non-zero elements that also contain zero elements from starting_param_subsets
    param_subsets = []
    for starting_param_subset_mask in starting_param_subset_masks
        push!(
            param_subsets,
            [
                subset_w_n_removed_params .* starting_param_subset_mask for
                subset_w_n_removed_params in subset_of_all_param_subsets_w_n_removed_params if
                sum(((subset_w_n_removed_params .* starting_param_subset_mask) .> 0)[1:end-n_non_params]) == n_removed_params
            ]...,
        )
    end
    param_subsets = [
        NamedTuple{Tuple(nt_param_choice_names)}(x) for
        x in unique(param_subsets) if sum((values(x) .> 0)[1:end-n_non_params]) == n_removed_params
    ]

    # Filter out parameter subsets that produce NaN loss values (e.g., when all substrate K = Inf)
    # add a new column to data to assign an integer to each source/figure from publication
    data = PKM2_data_for_fit
    data.fig_num = vcat(
        [
            i * ones(Int64, count(==(unique(data.source)[i]), data.source)) for
            i = 1:length(unique(data.source))
        ]...,
    )
    # convert DF to NamedTuple for better type stability / speed
    rate_data = Tables.columntable(
        data[.!isnan.(data.Rate), [:Rate, :PEP, :ADP, :Pyruvate, :ATP, :Phenylalanine, :F16BP, :fig_num]],
    )
    # make a vector containing indexes of points corresponding to each figure
    fig_point_indexes = [findall(data.fig_num .== i) for i in unique(data.fig_num)]
    # do the actual filtering
    param_subsets = [
        param_subset for param_subset in unique(param_subsets) if !isnan(
            loss_PKM2!(
                ones(length(param_names)),
                rate_data,
                fig_point_indexes;
                nt_param_choice = param_subset,
            ),
        )
    ]

    # Fit equation with param_subsets parameters to training data
    results_array = pmap(
        param_subset -> train_PKM2(data, param_names; n_iter = 20, nt_param_choice = param_subset),
        param_subsets,
    )

    # Save results of param combination, training loss, test loss, values of parameters
    column_names = [:n_removed_params, :error_counter, :param_subset, :train_loss, param_names...]
    data_for_table = permutedims(
        hcat(
            [
                [
                    n_removed_params
                    results_array[i][3]
                    string(param_subsets[i])
                    results_array[i][1]
                    param_rescaling(results_array[i][2])
                ] for i in eachindex(results_array)
            ]...,
        ),
    )
    df_results = DataFrame(data_for_table, column_names)
    CSV.write(
        "Cluster_results/$(Dates.format(now(),"mmddyy"))_PKM2_reverse_subset_selection_$(n_removed_params)_removed_params.csv",
        df_results,
    )

    # Filter subsets with NaN train or test loss values
    df_results = filter(row -> !isnan(row.train_loss), df_results)

    # Filter DataFrame df_result to only include rows that have train_loss lower than 1.1*minimum(train_loss)
    filter!(row -> row.train_loss < 1.1 * minimum(df_results.train_loss), df_results)

    # Sort DataFrame df_result by train_loss from low to high
    sort!(df_results, :train_loss)

    # Take the n_best_models param_subsets for this n_removed_params
    best_param_subsets_for_this_n_removed_params = eval.(Meta.parse.(df_results.param_subset))

    # Add the best param_subset for this n_removed_params to the list of best_param_subsets for later test loss calculation
    push!(best_param_subsets, best_param_subsets_for_this_n_removed_params...)

    # Update starting_param_subsets to be the best param_subsets for this n_removed_params
    global starting_param_subsets = values.(best_param_subsets_for_this_n_removed_params)
end





# Calculate test losses for each best model by training it on data missing one fig and testing on that fig

# n_repeats repeats of each fig and param_subset combination
n_repeats = 1

figs = unique(PKM2_data_for_fit.source)
list_figs = repeat(vcat([fill(fig, length(best_param_subsets)) for fig in figs]...), n_repeats)
list_param_subsets = repeat(vcat([best_param_subsets for fig in figs]...), n_repeats)
@assert length(list_figs) == length(list_param_subsets)

# Fit to training data missing one figure to rank subsets based on their training loss value
fig_jacknife_results_array = pmap(
    (fig, param_subset) ->
        fig_loocv_PKM2(fig, PKM2_data_for_fit, param_names; n_iter = 20, nt_param_choice = param_subset),
    list_figs,
    list_param_subsets,
)

# Save results of param combination, training loss, test loss, values of parameters
column_names =
    [:n_removed_params, :error_counter, :param_subset, :dropped_fig, :train_loss, :test_loss, param_names...]
data_for_table = permutedims(
    hcat(
        [
            [
                sum(values(list_param_subsets[i])[1:end-n_non_params] .> 0)
                fig_jacknife_results_array[i][5]
                list_param_subsets[i]
                fig_jacknife_results_array[i][1]
                fig_jacknife_results_array[i][2]
                fig_jacknife_results_array[i][3]
                param_rescaling(fig_jacknife_results_array[i][4])
            ] for i in eachindex(fig_jacknife_results_array)
        ]...,
    ),
)
df_results = DataFrame(data_for_table, column_names)
CSV.write(
    "Cluster_results/$(Dates.format(now(),"mmddyy"))_PKM2_reverse_subset_selection_train_test_loss_for_all_figs_$(n_repeats)_repeats.csv",
    df_results,
)



# Remove all the workers
for n in workers()
    rmprocs(n)
end
