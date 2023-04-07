#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=4G
#SBATCH --mail-user=msmoss@uwaterloo.ca
#SBATCH --mail-type=FAIL
#SBATCH --account=rrg-rgmelko-ab
#SBATCH --output=outputs/slurm-%j.out

module load StdEnv/2020 julia/1.8.1

julia run_MCMC_ising.jl \
   $L $Start_Chain