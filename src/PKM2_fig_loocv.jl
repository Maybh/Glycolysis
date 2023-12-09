# cd(@__DIR__)
# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()

# using Distributed, ClusterManagers, Dates

# # ENV["SLURM_JOB_CPUS_PER_NODE"] give number of cores in "40(x4),32,20" format
# # "40(x2),32,20" can be parsed to 132 using code below to get total number of allocated cores

# subs = Dict("x" => "*", "(" => "", ")" => "");
# np = sum(eval(Meta.parse(replace(ENV["SLURM_JOB_CPUS_PER_NODE"], r"x|\(|\)" => s -> subs[s]))))
# addprocs(SlurmManager(np); exeflags="--project")


#Uncomment code below to run on local machine or uncomment code above to run on cluster
#If run locally adjust number of workers to number of cores on your machine
using Distributed, Dates
addprocs(4; exeflags = "--project")

#= 
This code perform LOOCV to fit the rate_PKM2() equation to full data missing one of the figures.
This probably should be run on cluster as it will takes 1-3min to run one fit on local machine
and there are 20+ figures to fit to do full LOOCV
=#

@everywhere using CMAEvolutionStrategy, Optim, Random, Sobol, DataFrames, Statistics, CSV

@everywhere include("PKM2_fitting_functions.jl")

#Load and process data
PKM2_data_for_fit = CSV.read("PKM2_data/PKM2_data.csv", DataFrame)

#Add source column that uniquely identifies a figure from publication
PKM2_data_for_fit.source .= PKM2_data_for_fit.Article .* "_" .* PKM2_data_for_fit.Fig

#Filter the data to remove outliers
# filter!(
#     (row -> !(row.source == "Dombrauckas_Biochem_2005_Fig5A" && (row.PEP < 0.1e-3 || row.Rate < 0.2))),
#     PKM2_data_for_fit,
# )

#Train on data missing one figure and trest on that figure
n_repeats = 2 #n_repeats of fit without each figure
list_figs = vcat([fill(fig, n_repeats) for fig in unique(PKM2_data_for_fit.source)]...)

fig_loocv_results_array =
    @time pmap(x -> fig_loocv_PKM2(x, PKM2_data_for_fit, param_names; n_iter = 20), list_figs)
column_names = [:dropped_fig, :train_fit, :test_fit, param_names...]
data_for_table = permutedims(
    hcat(
        [
            [
                fig_loocv_results_array[i][1]
                fig_loocv_results_array[i][2]
                fig_loocv_results_array[i][3]
                param_rescaling(fig_loocv_results_array[i][4])
            ] for i in eachindex(list_figs)
        ]...,
    ),
)

CSV.write(
    "Cluster_results/$(Dates.format(now(),"mmddyy"))_PKM2_fig_loocv_$(n_repeats)_fig_repeats.csv",
    DataFrame(data_for_table, column_names),
)

#Train on all data and test on one figure
bootstrap_param_data = CSV.read(
    "Results/102523_PKM2_fitting_8_not_rescaled_params_full_model.csv",
    DataFrame,
)
kinetic_params = Vector(
    eachrow(
        bootstrap_param_data[:, param_names][
            bootstrap_param_data.fitness.==minimum(bootstrap_param_data.fitness),
            :,
        ],
    )[1],
)
fig_loocv_results_array_train_all_fig =
    @time pmap(x -> test_PKM2(PKM2_data_for_fit[PKM2_data_for_fit.source.==x, :], kinetic_params), list_figs)

column_names = [:dropped_fig, :test_fit]
data_for_table = permutedims(
    hcat(
        [
            [
                list_figs[i]
                fig_loocv_results_array_train_all_fig[i]
            ] for i in eachindex(list_figs)
        ]...,
    ),
)
CSV.write(
    "Cluster_results/$(Dates.format(now(),"mmddyy"))_PKM2_fig_loocv_full_train_w_all_data.csv",
    DataFrame(data_for_table, column_names),
)

# Remove all the workers
for n in workers()
    rmprocs(n)
end
