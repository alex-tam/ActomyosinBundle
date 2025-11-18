# Control function for repeated actomyosin network simulations
# Alex Tam, 11/07/2022

# Import packages
using Revise # Reduce need to restart Julia      
using Parameters # Tools for storing parameters in data structures
using LinearAlgebra # Matrix and vector operations
using StatsBase # Sampling without replacement
using Plots # Plotting library
using Plots.PlotMeasures # Enable re-sizing plot margins
using LaTeXStrings # Display LaTeX output
using ForwardDiff # Automatic differentiation tools
using Optim # Optimisation routines
using LineSearches # Line searches for optimisation
using Printf # C-style printing macros
using Distributions # Probability distributions
using DelimitedFiles # Read and write delimited data
using Loess # LOESS regression

# Include code from external files
include("model_parameters.jl")
include("actomyosin_network.jl")
include("State.jl")
include("initial_condition.jl")
include("Actin_Filament.jl")
include("periodic.jl")
include("Cross_Link.jl")
include("Myosin_Motor.jl")
include("intersection_search.jl")
include("thermal.jl")
include("network_statistics.jl")
include("network_force.jl")
include("draw_network.jl")
include("turnover.jl")
include("optimise_network.jl")
include("energy.jl")
include("gradient.jl")
include("draw_statistics.jl")

function run_simulations()
    ##### Establish plots #####
    gr() # Load GR plotting back-end
    default(titlefont = (14, "Computer Modern"), guidefont = (18, "Computer Modern"), tickfont = (14, "Computer Modern")) # Plot settings
    ##### Paremeters to simulate #####
    params = [0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20] # Vector of parameter values to simulate
    ##### Loop over parameter values #####
    for par in params 
        ##### Set parameters #####
        parN = Numerical_Parameters(mechanics = true) # Initialise struct of numerical parameters
        parA = Actin_Properties(k_uniform = 0.0, k_branching = par, k_bias = 0.0, k_t = 0.0) # Initialise struct of actin filament properties
        parM = Myosin_Properties() # Initialise struct of myosin motor properties
        writedlm("t.csv", (0:parN.nT-1).*parN.dt) # Output vector of discretised time
        ##### Pre-allocate time-averaged quantities (all trials) ####
        stress = Vector{Float64}() # Pre-allocate time-averaged mean normal stress
        Ia = Vector{Float64}() # Pre-allocate time-averaged aggregation statistic
        Ib = Vector{Float64}() # Pre-allocate time-averaged bundle statistic
        Ip = Vector{Float64}() # Pre-allocate time-averaged parallel bundle statistic
        mp = Vector{Float64}() # Pre-allocate time-averaged mean relative motor position
        K = Vector{Float64}() # Pre-allcoate time-averaged mean integrated filament curvature
        θ = Vector{Float64}() # Pre-allcoate time-averaged mean motor angle
        I2 = Vector{Float64}() # Pre-allcoate time-averaged mean motor angle
        for j = 1:parN.nTrials # Repeat trials
            ##### Run simulations ####
            NetStats = actomyosin_network(parN, parA, parM, par, j, false)
            ##### Store time-averaged quantities (all trials) #####
            push!(stress, NetStats.σ) # Time-averaged mean normal stress
            push!(Ia, NetStats.Ia) # Time-averaged aggregation statistic
            push!(Ib, NetStats.Ib) # Time-averaged bundle statistic
            push!(Ip, NetStats.Ip) # Time-averaged parallel bundle statistic
            push!(mp, NetStats.mp) # Time-averaged mean relative motor position
            push!(K, NetStats.K) # Time-averaged mean integrated filament curvature
            push!(θ, NetStats.θ) # Time-averaged mean motor angle
            push!(I2, NetStats.I2) # Time-averaged mean motor angle
            ##### Write trial data to files #####
            # Time-dependent (single trial)
            writedlm("mean_normal_stress-par-$par-trial-$j.csv", NetStats.Mean_Normal_Stress)
            writedlm("aggregation_statistic-par-$par-trial-$j.csv", NetStats.Aggregation_Statistic)
            writedlm("bundle_statistic-par-$par-trial-$j.csv", NetStats.Bundle_Statistic)
            writedlm("parallel_statistic-par-$par-trial-$j.csv", NetStats.Parallel_Statistic)
            writedlm("mean_motor_position-par-$par-trial-$j.csv", NetStats.Motor_Position)
            writedlm("mean_integrated_curvature-par-$par-trial-$j.csv", NetStats.Curvature)
            writedlm("mean_motor_angle-par-$par-trial-$j.csv", NetStats.Motor_Angle)
            writedlm("mean_2f_index-par-$par-trial-$j.csv", NetStats.TwoFilamentIndex)
            # Time-averaged (all trials)
            writedlm("mean_normal_stress_time-averaged-par-$par.csv", stress)
            writedlm("aggregation_statistic_time-averaged-par-$par.csv", Ia)
            writedlm("bundle_statistic_time-averaged-par-$par.csv", Ib)
            writedlm("parallel_statistic_time-averaged-par-$par.csv", Ip)
            writedlm("mean_motor_position_time-averaged-par-$par.csv", mp)
            writedlm("mean_integrated_curvature_time-averaged-par-$par.csv", K)
            writedlm("mean_motor_angle_time-averaged-par-$par.csv", θ)
            writedlm("mean_2f_index_time-averaged-par-$par.csv", I2)
        end
    end
end

@time run_simulations()
