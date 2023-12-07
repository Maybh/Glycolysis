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
This code fits the rate_PKM2() equation to 
either full data 
or bootstrap data consisting of random figures taking with replacement
or bootstrap data consisting of random point taking with replacement
Uncomment the code in the bottom to choose between options above
=#

@everywhere using CMAEvolutionStrategy, Optim, Random, Sobol, DataFrames, Statistics, CSV, StatsBase

#Load and process data
PKM2_data_for_fit = CSV.read("PKM2_data/PKM2_data.csv", DataFrame)

#Add source column that uniquely identifies a figure from publication
PKM2_data_for_fit.source .= PKM2_data_for_fit.Article .* "_" .* PKM2_data_for_fit.Fig

#Filter the data to remove outliers
# filter!(
#     (row -> !(row.source == "Dombrauckas_Biochem_2005_Fig5A" && (row.PEP < 0.1e-3 || row.Rate < 0.2))),
#     PKM2_data_for_fit,
# )

@everywhere include("PKM2_fitting_functions.jl")

# Uncomment below to fit to full data
data = PKM2_data_for_fit
number_bootstrap_datasets = 8
bootstrap_datasets = [data for i = 1:number_bootstrap_datasets]

# # Uncomment below to fit to bootstrap datasets from points
# data = PKM2_data_for_fit
# number_bootstrap_datasets = 10_000
# bootstrap_datasets = [
#     data[sample(axes(data, 1), nrow(data); replace = true, ordered = true), :] for
#     i = 1:number_bootstrap_datasets
# ]

# # Uncomment below to fit to bootstrap datasets from figures
# data = PKM2_data_for_fit
# list_figs = unique(PKM2_data_for_fit.source)
# number_bootstrap_datasets = 10_000
# number_repeats = 5
# bootstrap_fig_list = vcat([fill(sample(list_figs, length(list_figs), replace=true), number_repeats) for i in 1:number_bootstrap_datasets]...)
# bootstrap_datasets = [vcat([data[data.source .== fig, :] for fig in bootstrap_fig]...) for bootstrap_fig in bootstrap_fig_list]

# Fit bootstrapped datasets
res_array = @time pmap(x -> train_PKM2(x, param_names; n_iter = 20), bootstrap_datasets)

# print(res_array[1][3])

# Save bootstrapping results
data_for_table = permutedims(
    hcat([[
        # 1 + (i - 1) ÷ number_repeats
        # [bootstrap_fig_list[i]]
        #uncomment two lines above for fig bootstrap
        res_array[i][1]
        # param_rescaling(res_array[i][2])
        #uncomment 1 line above to save non rescaled params for testing
        #uncomment 1 line below to save non rescaled params for testing
        res_array[i][2]
    ] for i in eachindex(res_array)]...),
)
column_names = [
    # :bootrap_number,
    # :bootrap_fig_list, 
    #uncomment two lines above for fig bootstrap
    :fitness,
    param_names...,
]

Data_for_CSV = DataFrame(data_for_table, column_names)
CSV.write(
    "Results/$(Dates.format(now(),"mmddyy"))_PKM2_fitting_$(number_bootstrap_datasets)_not_rescaled_params_full_model.csv",
    Data_for_CSV,
)

println(minimum(Data_for_CSV.fitness))

# Remove all the workers
for n in workers()
    rmprocs(n)
end


