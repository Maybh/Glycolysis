#!/bin/bash
# Job name:
#SBATCH --job-name=julia_bboptim
#
# Account:
#SBATCH --account=fc_moorjani
#
# Partition:
#SBATCH --partition=savio3
#
# Request one node:
#SBATCH --nodes=24
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=10:00:00

# Load software
module purge
module load julia/1.9.2

# Run Julia script
julia PKM2_reverse_model_selection.jl
