# #!/bin/bash

# FSS Dataset

joblist=$(sq -h --format="%j")
for L in 8 12 16 20 24 32 
do
    X="q1b|$L|chains1-5"
    sbatch -J "$X" --export="L=$L,Start_Chain=1" submit_MCMC_ising.sh

    X="q1b|$L|chains6-10"
    sbatch -J "$X" --export="L=$L,Start_Chain=6" submit_MCMC_ising.sh

    X="q1b|$L|chains10-15"
    sbatch -J "$X" --export="L=$L,Start_Chain=11" submit_MCMC_ising.sh

    X="q1b|$L|chains16-20"
    sbatch -J "$X" --export="L=$L,Start_Chain=16" submit_MCMC_ising.sh

    sleep 0.5s
done
