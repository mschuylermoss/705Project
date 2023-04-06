using Random
using DataFrames
using CSV

function build_lattice_df(Lx,Ly)
    
    N = Lx*Ly
    interactions = DataFrame(["$spin" => [] for spin in 1:N])    
    for spin in 1:N
        up = mod(spin-1-Lx,N)
        down = mod(spin-1+Lx,N)
        right = mod(spin-1+1,N)
        left = mod(spin-1-1,N)
        push!(interactions[!,"$spin"],[up+1,down+1,left+1,right+1])
    end
    
    return interactions
end

function Calculate_Ising_EnergyPerSpin(config,interactions,J)
    
    N = length(config)
    energy = 0
    
    for spin in (1:N)
        neighbors = interactions[!,"$spin"][1,1]
        energy_contributions = 0.5*-1*J*config[spin].*config[neighbors]
        energy += sum(energy_contributions)
    end
    
    return energy/N
    
end

function build_exp_beta_deltaE_df(beta)
    
    possible_sum_neighbors = collect(-4:2:4)
    spins = [1,-1]
    exp_beta_deltaEs = DataFrame(["$spin,$poss" => [] for poss in possible_sum_neighbors for spin in spins])    
    for sum_neighbors in possible_sum_neighbors
        up_exp_val = exp(-2*beta*sum_neighbors)
        down_exp_val = exp(2*beta*sum_neighbors)
        push!(exp_beta_deltaEs[!,"1,$sum_neighbors"],up_exp_val)
        push!(exp_beta_deltaEs[!,"-1,$sum_neighbors"],down_exp_val)
    end
    
    return exp_beta_deltaEs
end

function Metropolis_update_function(config,interactions,beta,J,exp_beta_deltaE_vals)
    
    rand_spin_flip_index = rand(1:length(config),1)[1]
    new_config = copy(config)

    spin = Int(config[rand_spin_flip_index])
    neighbors = interactions[!,"$rand_spin_flip_index"][1,1]
    sum_neighbors = Int(sum(config[neighbors]))
    exp_beta_deltaE = exp_beta_deltaE_vals[!,"$spin,$sum_neighbors"][1]
    rand_value = rand(Float64,1)[1]
    
    if exp_beta_deltaE > rand_value
        new_config[rand_spin_flip_index] = -1*config[rand_spin_flip_index]
    end
    
    return new_config
end

function MonteCarlo_Ising(beta,J,Lx,Ly,warmup_steps,steps,build_lattice_function,build_update_probs,update_function,energy_function)
    
    # system information
    N = Lx*Ly
    interactions = build_lattice_function(Lx,Ly)
    exp_beta_deltaEs = build_exp_beta_deltaE_df(beta)
    initial_config = rand((-1,1),N)
    
    # initialize dataframes
    fulldf = DataFrame(:Es => Float64[],:E2s => Float64[],:Ms => Int[],:M2s => Int[])
    maindf = DataFrame(:avenergy => Float64[],:avenergy2 => Float64[],:avM => Float64[],:avM2 => Float64[])
    
    # warmup steps - track Es and Ms but not in averages
    config = initial_config
    for step in (1:warmup_steps)
        new_config = update_function(config,interactions,beta,J,exp_beta_deltaEs)
        new_energy = energy_function(new_config,interactions,J)
        M = sum(new_config)
        push!(fulldf,(new_energy,new_energy^2,M,M^2))
        config = copy(new_config)
    end
    
    # begin MC steps and observable tracking
    sumE = 0
    sumE2 = 0
    sumM = 0
    sumM2 = 0
    println("Warmup steps over!")
    for step in (1:steps)
        new_config = update_function(config,interactions,beta,J,exp_beta_deltaEs)
        new_energy = energy_function(new_config,interactions,J)
        M = sum(new_config)
        push!(fulldf,(new_energy,new_energy^2,M,M^2))
        absM = abs(M)
        sumE = sumE+new_energy
        sumE2 = sumE2+new_energy^2
        sumM = sumM+absM
        sumM2 = sumM2+absM^2
        push!(maindf, (sumE/step,sumE2/step,sumM/step,sumM2/step))
        config = copy(new_config)
    end
        
    return fulldf, maindf
end