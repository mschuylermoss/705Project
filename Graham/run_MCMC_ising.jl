using Random
using DataFrames
using CSV
using ArgParse
include("MCMC_ising.jl")

s = ArgParseSettings()
@add_arg_table! s begin
    "L"
        help = "The length of the square lattice (PBC dimension)"
        required = true
        arg_type = Int
    "Start_Chain"
        help = "The first chain to run"
        required = true
        arg_type = Int
end
parsed_args = parse_args(ARGS, s)

L = parsed_args["L"]
Start_Chain = parsed_args["Start_Chain"]
Ts = 1.0:0.05:4.0
J = 1
N_chains = 5
n_steps = 1000000
warmup = Int(n_steps/10)

for T in Ts
    beta = 1/T
    T = round(T,digits=2)
    for chain in (Start_Chain):(Start_Chain+N_chains-1)
        inner_steps = min((L*L)/2,200)
        avgs = MonteCarlo_Ising(beta,J,L,L,warmup,n_steps,inner_steps,build_lattice_df,build_deltaE_df,Metropolis_update_function,Calculate_Ising_Energy);
        CSV.write("/scratch/msmoss/705_Project/q1b_newcode/L_$(L)/T_$(T)_chain$(chain)_avgvals.csv",  avgs, writeheader=true)
    end
end
