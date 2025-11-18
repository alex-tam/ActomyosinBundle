# Data structures for model parameters
# Alex Tam, 13/10/2020

# Parameters
"Global numerical parameters"
@with_kw struct Numerical_Parameters
    nA::Int = 150 # [-] Number of actin filaments
    nM::Int = 30 # [-] Number of myosin motors
    nT::Int = 6001 # [-] Number of time steps
    nTrials::Int = 10 # [-] Number of trials for each parameter value
    plot_interval::Int = 1200 # [-] Time steps per network plot
    dt::Float64 = 0.05 # [s] Time step size
    lxx::Float64 = 5.0 # [μm] Reference domain width (x)
    lyy::Float64 = 5.0 # [μm] Reference domain width (y)
    xTol::Float64 = 1e-8 # [-] DOF tolerance for optimisation
    fTol::Float64 = 1e-8 # [-] Objective function tolerance for optimisation
    gTol::Float64 = 1e-8 # [-] Gradient tolerance for optimisation
    mechanics::Bool = true
end

"Actin filament properties"
@with_kw struct Actin_Properties
    nSeg::Int = 5 # [-] Number of segments
    LSeg::Float64 = 0.2 # [μm] Equilibrium segment length
    k::Float64 = 1000 # [pN/μm] Spring constant
    lambda_a::Float64 = 0.05 # [pN/(μm^2)*s] Actin-background drag coefficient
    kappa::Float64 = 0.073 # [pN*μm^2] Flexural rigidity
    lambda_pf::Float64 = 30.0 # [pN/μm*s] Protein friction drag coefficient
    k_uniform::Float64 = 0.04 # [/filament/s] Uniform turnover rate
    k_bias::Float64 = 0.04 # [/filament/s] Biased turnover rate
    k_branching::Float64 = 0.04 # [/filament/s] Branching turnover rate
    k_t::Float64 = 0.0 # [/filament/s] Treadmilling rate
    kb::Float64 = 1.380649e-5 # [μm*pN/K] Boltzmann constant
    T::Float64 = 298.15 # [K] Temperature
    σ_p::Float64 = 0.05 # [-] Standard deviation for translations in biased turnover
    σ_a::Float64 = π/20 # [-] Standard deviation for rotations in biased turnover
    PRev::Float64 = 0.5 # [-] Probability of reversed polarity (biased turnover)
    θb::Float64 = 70 # [degrees] Branching angle
end

"Myosin motor properties"
@with_kw struct Myosin_Properties
    k::Float64 = 1000 # [pN/μm] Spring constant
    Fs::Float64 = 5 # [pN] Stall force
    Vm::Float64 = 0.5 # [μm/s] Load-free velocity
    k_off::Float64 = 0.35 # [/motor/s] Reference off-rate in the absence of force
    F_ref::Float64 = 12.6 # [pN] Reference force for computing force-dependent off-rate
end