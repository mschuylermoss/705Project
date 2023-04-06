# #!/bin/bash

# FSS Dataset

joblist=$(sq -h --format="%j")
for L in 8 12 16 20 24 32 #$(seq 100000 50000 200000)
do
    X="q1b|$L"
    sbatch -J "$X" --export="L=$L" submit_MCMC_ising.sh
    sleep 0.5s
done
