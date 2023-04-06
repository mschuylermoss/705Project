using Random
using DataFrames
using CSV
using ArgParse
include("MCMC_ising.jl")

s = ArgParseSettings()
@add_arg_table! s["groundstate"] begin
    "L"
        help = "The length of the square lattice (PBC dimension)"
        required = true
        arg_type = Int
end
parsed_args = parse_args(ARGS, s)

L = parsed_args["L"]
Ts = 1.0:0.05:4.0
J = 1
N_chains = 10
n_steps = 1000000
warmup = Int(n_steps/10)

for T in Ts
    beta = 1/T
    T = round(T,digits=2)
    for chain in 0:N_chains-1
        full,avgs = MonteCarlo_Ising(beta,J,L,L,warmup,n_steps,build_lattice_df,build_exp_beta_deltaE_df,Metropolis_update_function,Calculate_Ising_EnergyPerSpin);
        CSV.write("./data/q1b/L_$(L)/allvals_T_$(T)_chain$(chain).csv",  full, writeheader=true)
        CSV.write("./data/q1b/L_$(L)/avgvals_T_$(T)_chain$(chain).csv",  avgs, writeheader=true)
    end
end
