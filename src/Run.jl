cd(@__DIR__)
using Pkg
# cd("Glycolysis")
Pkg.activate(".")
# Pkg.activate("/global/home/users/maybenhamo/.julia/environments/v1.9")
dir_path = "/global/home/users/maybenhamo/Glycolysis"


@info "starting"

# using Distributed, ClusterManagers, Dates, JSON

# #Increase time before worker terminates without hearing from master process. Useful for large number of cores.
# ENV["JULIA_WORKER_TIMEOUT"] = 600.0

# # ENV["SLURM_JOB_CPUS_PER_NODE"] give number of cores in "40(x4),32,20" format
# # "40(x2),32,20" can be parsed to 132 using code below to get total number of allocated cores

# subs = Dict("x" => "*", "(" => "", ")" => "");
# np = sum(eval(Meta.parse(replace(ENV["SLURM_JOB_CPUS_PER_NODE"], r"x|\(|\)" => s -> subs[s]))))
# @info "running with #$np workers"
# addprocs(SlurmManager(np); exeflags = "--project")

using Distributed, ClusterManagers, Dates

include("Utils.jl")

config_path = joinpath(dir_path, "resource/config.json")
config = read_config(config_path)

@info "config: $config"

if lowercase(config["mode"]) == "cluster"
    initial_cluster()
elseif lowercase(config["mode"]) == "local"
    initial_local_distribution(n_process=16)
else
    throw(ArgumentError("mode is missing in configuration or invalid, it should be cluster or local]"))
end

@everywhere using Logging
@everywhere include("Utils.jl")
@everywhere include("VarSelection.jl")

global date_time = Dates.format(now(), "mmddyy HH:MM:SS")


function Run(config, data_path, dir_path)
    starting_time = now()
    formatted_time = Dates.format(starting_time, "HH:MM:SS")
    prefix = joinpath(dir_path, "logs/julia-$formatted_time")

    try
        # Load data
        global data = load_prepare_data(data_path)
        variable_selection(config, data, dir_path)
    catch e
        @error "An error occurred in variable_selection" exception=e
        throw(e)
    finally 
        rmprocs(workers())
    end
    
    # redirect_stdio(stdout=prefix * ".log", stderr=prefix * ".log", stdin=devnull) do
    #     try
    #         # Load data
    #         global data = load_prepare_data(data_path)
    #         variable_selection(config, data)
    #     catch e
    #         @error "An error occurred in variable_selection" exception=e
    #         throw(e)
    #     finally 
    #         rmprocs(workers())
    #     end
    # end
    
end

Run(config, joinpath(dir_path, "PKM2_data/PKM2_data.csv"), dir_path)