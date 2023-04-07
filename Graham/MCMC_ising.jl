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

function Calculate_Ising_Energy(config,interactions,J)
    
    N = length(config)
    energy = 0
    
    for spin in (1:N)
        neighbors = interactions[!,"$spin"][1,1]
        energy_contributions = 0.5*-1*J*config[spin].*config[neighbors]
        energy += sum(energy_contributions)
    end
    
    return energy
    
end

function build_deltaE_df(beta)
    
    possible_sum_neighbors = collect(-4:2:4)
    spins = [1,-1]
    deltaEs = DataFrame(["$spin,$poss" => [] for poss in possible_sum_neighbors for spin in spins])    
    for sum_neighbors in possible_sum_neighbors
        spinup_deltaE = 2*sum_neighbors
        spindown_deltaE = -2*sum_neighbors 
        push!(deltaEs[!,"1,$sum_neighbors"],spinup_deltaE)
        push!(deltaEs[!,"-1,$sum_neighbors"],spindown_deltaE)
    end
    
    return deltaEs
end

function Metropolis_update_function(config,interactions,beta,J,deltaE_vals)
    
    rand_spin_flip_index = rand(1:length(config),1)[1]
    new_config = copy(config)

    spin = Int(config[rand_spin_flip_index])
    neighbors = interactions[!,"$rand_spin_flip_index"][1,1]
    sum_neighbors = Int(sum(config[neighbors]))
    deltaE = deltaE_vals[!,"$spin,$sum_neighbors"][1]
    exp_beta_deltaE = exp(-1*beta*deltaE)
    rand_value = rand(Float64,1)[1]
    
    if exp_beta_deltaE > rand_value
        new_config[rand_spin_flip_index] = -1*config[rand_spin_flip_index]
    else
        deltaE = 0.
    end
    
    return new_config, deltaE
end

function MonteCarlo_Ising(beta,J,Lx,Ly,warmup_steps,steps,steps_per_step,build_lattice_function,build_update_probs,update_function,energy_function)
    
    # system information
    N = Lx*Ly
    interactions = build_lattice_function(Lx,Ly)
    deltaEs = build_deltaE_df(beta)
    initial_config = rand((-1,1),N)
    
    # initialize dataframes
    maindf = DataFrame(:avenergy => Float64[],:avenergy2 => Float64[],:avM => Float64[],:avM2 => Float64[])
    
    # warmup steps 
    config = initial_config
    energy = energy_function(config,interactions,J)
    for step in (1:warmup_steps)
        for inner_steps in (1:steps_per_step)
            new_config, deltaE = update_function(config,interactions,beta,J,deltaEs)
            config = copy(new_config)
            new_energy = energy+deltaE
            energy = new_energy
        end
    end
    
    # begin MC steps and observable tracking
    sumE = 0
    sumE2 = 0
    sumM = 0
    sumM2 = 0
    println("Warmup steps over!")
    for step in (1:steps)
        for inner_step in (1:steps_per_step)
            new_config,deltaE = update_function(config,interactions,beta,J,deltaEs)
            config = copy(new_config)
            new_energy = energy+deltaE
            energy = new_energy
        end
        M = sum(config)
        absM = abs(M)
        sumE = sumE+energy
        sumE2 = sumE2+energy^2
        sumM = sumM+absM
        sumM2 = sumM2+absM^2
        push!(maindf, (sumE/step,sumE2/step,sumM/step,sumM2/step))
    end
        
    return maindf
end