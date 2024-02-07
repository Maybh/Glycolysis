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
# Write logs:
#SBATCH --output=/global/home/users/maybenhamo/Glycolysis/logs/slurm-%j.out # STDOUT
#SBATCH --error=/global/home/users/maybenhamo/Glycolysis/logs/slurm-%j.err
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=24:00:00

# Load software
module purge
module load julia/1.9.2

# Run Julia script
julia src/Run.jl --project=/global/home/users/maybenhamo/.julia/environments/v1.9