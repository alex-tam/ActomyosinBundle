# Simulate a 2D actomyosin network
# Alex Tam, 11/07/2022

"Perform one actomyosin network simulation"
function actomyosin_network(parN, parA, parM, par, trial, animate::Bool)
    # Specify physical patch dimensions [um]
    Lxx::Float64 = 5.0 # Width (x)
    Lxy::Float64 = 0 # Zero for rectangular patch
    Lyx::Float64 = 0 # Zero for rectangular patch
    Lyy::Float64 = 5.0 # Height (y)
    # Pre-allocate output vectors (time-varying)
    Force = [[0.0, 0.0, 0.0, 0.0] for idx in 1:parN.nT] # [pN] Force vectors (Fxx, Fxy, Fyx, Fyy) for all time steps
    Aggregation_Statistic = [0.0 for idx in 1:parN.nT] # [-] Aggregation index
    Bundle_Statistic = [0.0 for idx in 1:parN.nT] # [-] Bundle formation index
    Parallel_Statistic = [0.0 for idx in 1:parN.nT] # [-] Parallel bundle index
    Bundle_Direction = [0.0 for idx in 1:parN.nT] # [rad] Principal direction of fibres in bundled structure
    Motor_Position = [0.0 for idx in 1:parN.nT] # [-] Mean relative motor position
    Curvature = [0.0 for idx in 1:parN.nT] # [-] Mean integrated filament curvature
    Motor_Angle = [0.0 for idx in 1:parN.nT] # [-] Mean motor angle
    TwoFilamentIndex = [0.0 for idx in 1:parN.nT] # [-] Mean two-filament index
    # Generate initial conditions
    state = State{Float64}(Vector{Vector{Vector}}(), Vector{Vector}()) # Initialise empty State struct
    mm = Vector{Myosin_Motor}() # Pre-allocate empty myosin motors
    af, state = actin_ic(state, parN, parA, Lxx, Lyy) # Initialise actin filaments
    xl = intersection_search(state, af, mm) # Initialise cross-links
    mm, xl, state = myosin_ic(state, mm, parN, xl, Lxx, Lxy, Lyx, Lyy) # Initialise myosin motors
    state_old = state # Store initial state to compute initial force vector
    random = thermal(af, state) # Random variables for thermal motion
    # Time-stepping
    if parN.mechanics == true # If mechanical model is included
        if animate == true # Produce animation
            animation = @animate for i = 1:parN.nT
                # Write final DOF vector to file
                if i == parN.nT
                    dof = build_dof(state)
                    writedlm("dof-par-$par-trial-$trial-step-$i.csv", dof)
                end
                # Compute network statistics
                Force[i] = network_force(state, state_old, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) # Force vector
                Aggregation_Statistic[i], Bundle_Statistic[i], Parallel_Statistic[i], bundled_directions = bundle_statistic(state, af) # Compute bundle structure statistics
                Bundle_Direction[i] = mean(bundled_directions) # Compute principle bundle direction
                Motor_Position[i] = mean_relative_mp(state)
                Curvature[i] = curvature(af, state, Lxx, Lxy, Lyx, Lyy)
                Motor_Angle[i], TwoFilamentIndex[i] = two_filament_index(mm, state, parA, Lxx, Lxy, Lyx, Lyy)
                # Draw network
                if mod(i-1, parN.plot_interval) == 0 
                    # draw_bundle_histogram(bundled_directions, Force[i], Lxx, Lyy, par, trial, i) # Draw and save bundle histogram
                    dof = build_dof(state)
                    writedlm("dof-par-$par-trial-$trial-step-$i.csv", dof)
                end
                draw_network(state, af, xl, mm, parN, parA, Lxx, Lxy, Lyx, Lyy, Force[i], false)
                if mod(i-1, parN.plot_interval) == 0 
                    savefig("actomyosin-par-$par-trial-$trial-step-$i.pdf") # Save network image
                end
                # Simulate one time step
                if i != parN.nT # Ensure correct looping sequence
                    af, mm, state = actin_turnover(state, af, mm, parN, parA, Lxx, Lyy) # Implement filament turnover
                    af, mm, state = actin_treadmilling(state, af, mm, parN, parA) # Implement filament treadmilling
                    af = segment_translations(state, af) # Update periodic BC information
                    xl = intersection_search(state, af, mm) # Update cross-links
                    mm, xl, state = myosin_turnover(state, xl, mm, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) # Myosin binding/unbinding
                    random = thermal(af, state) # Random variables for thermal motion
                    # Compute solution and store data
                    state_old = state; # Store current state for energy functional
                    new_dof = optimise_network(state_old, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) # Compute network solution
                    state = build_state(new_dof, af, mm) # Construct State from vector
                end
            end
            # Post-processing
            gif(animation, "actomyosin-par-$par-trial-$trial.gif", fps = 10)
            NetStats = network_statistics(parN, Force, Aggregation_Statistic, Bundle_Statistic, Parallel_Statistic, Bundle_Direction, Motor_Position, Curvature, Motor_Angle, TwoFilamentIndex, Lxx, Lyy)
            draw_statistics(NetStats, parN, parN.nT, par, trial)
        else # No animation
            for i = 1:parN.nT
                # Write final DOF vector to file
                if i == parN.nT
                    dof = build_dof(state)
                    writedlm("dof-par-$par-trial-$trial-step-$i.csv", dof)
                end
                # Compute network statistics
                Force[i] = network_force(state, state_old, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) # Force vector
                Aggregation_Statistic[i], Bundle_Statistic[i], Parallel_Statistic[i], bundled_directions = bundle_statistic(state, af) # Compute bundle structure statistics
                Bundle_Direction[i] = mean(bundled_directions) # Compute principle bundle direction
                Motor_Position[i] = mean_relative_mp(state)
                Curvature[i] = curvature(af, state, Lxx, Lxy, Lyx, Lyy)
                Motor_Angle[i], TwoFilamentIndex[i] = two_filament_index(mm, state, parA, Lxx, Lxy, Lyx, Lyy)
                # Draw network
                if mod(i-1, parN.plot_interval) == 0
                    draw_network(state, af, xl, mm, parN, parA, Lxx, Lxy, Lyx, Lyy, Force[i], false)
                    savefig("actomyosin-par-$par-trial-$trial-step-$i.pdf") # Save network image
                    dof = build_dof(state)
                    writedlm("dof-par-$par-trial-$trial-step-$i.csv", dof)
                    # draw_bundle_histogram(bundled_directions, Force[i], Lxx, Lyy, par, trial, i) # Draw and save bundle histogram
                end
                # Simulate one time step
                if i != parN.nT # Ensure correct looping sequence
                    af, mm, state = actin_biased_turnover(state, af, mm, parN, parA, Lxx, Lyy) # Implement biased filament turnover
                    af, mm, state = actin_branching_turnover(state, af, mm, parN, parA, Lxx, Lyy) # Implement filament branching
                    af, mm, state = actin_uniform_turnover(state, af, mm, parN, parA, Lxx, Lyy) # Implement uniform filament turnover
                    af, mm, state = actin_treadmilling(state, af, mm, parN, parA) # Implement filament treadmilling
                    af = segment_translations(state, af) # Update periodic BC information
                    xl = intersection_search(state, af, mm) # Update cross-links
                    mm, xl, state = myosin_turnover(state, xl, mm, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) # Myosin binding/unbinding
                    random = thermal(af, state) # Random variables for thermal motion
                    # Compute solution and store data
                    state_old = state; # Store current state for energy functional
                    new_dof = optimise_network(state_old, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) # Compute network solution
                    state = build_state(new_dof, af, mm) # Construct State from vector
                end
            end
            # Post-processing
            NetStats = network_statistics(parN, Force, Aggregation_Statistic, Bundle_Statistic, Parallel_Statistic, Bundle_Direction, Motor_Position, Curvature, Motor_Angle, TwoFilamentIndex, Lxx, Lyy)
            draw_statistics(NetStats, parN, parN.nT, par, trial)
        end
    else # No mechanical model
        if animate == true # Produce animation
            animation = @animate for i = 1:parN.nT
                # Write final DOF vector to file
                if i == parN.nT
                    dof = build_dof(state)
                    writedlm("dof-par-$par-trial-$trial-step-$i.csv", dof)
                end
                # Compute network statistics
                Force[i] = network_force(state, state_old, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) # Force vector
                Aggregation_Statistic[i], Bundle_Statistic[i], Parallel_Statistic[i], bundled_directions = bundle_statistic(state, af) # Compute bundle structure statistics
                Bundle_Direction[i] = mean(bundled_directions) # Compute principle bundle direction
                Motor_Position[i] = mean_relative_mp(state)
                Curvature[i] = curvature(af, state, Lxx, Lxy, Lyx, Lyy)
                Motor_Angle[i], TwoFilamentIndex[i] = two_filament_index(mm, state, parA, Lxx, Lxy, Lyx, Lyy)
                # Draw network
                if mod(i-1, parN.plot_interval) == 0 
                    # draw_bundle_histogram(bundled_directions, Force[i], Lxx, Lyy, par, trial, i) # Draw and save bundle histogram
                    dof = build_dof(state)
                    writedlm("dof-par-$par-trial-$trial-step-$i.csv", dof)
                end
                draw_network(state, af, xl, mm, parN, parA, Lxx, Lxy, Lyx, Lyy, Force[i], false)
                if mod(i-1, parN.plot_interval) == 0 
                    savefig("actomyosin-par-$par-trial-$trial-step-$i.pdf") # Save network image
                end
                # Simulate one time step
                if i != parN.nT # Ensure correct looping sequence
                    af, mm, state = actin_turnover(state, af, mm, parN, parA, Lxx, Lyy) # Implement filament turnover
                    af, mm, state = actin_treadmilling(state, af, mm, parN, parA) # Implement filament treadmilling
                    af = segment_translations(state, af) # Update periodic BC information
                    xl = intersection_search(state, af, mm) # Update cross-links
                    mm, xl, state = myosin_turnover(state, xl, mm, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) # Myosin binding/unbinding
                end
            end
            # Post-processing
            gif(animation, "actomyosin-par-$par-trial-$trial.gif", fps = 10)
            NetStats = network_statistics(parN, Force, Aggregation_Statistic, Bundle_Statistic, Parallel_Statistic, Bundle_Direction, Motor_Position, Curvature, Motor_Angle, TwoFilamentIndex, Lxx, Lyy)
            draw_statistics(NetStats, parN, parN.nT, par, trial)
        else # No animation
            for i = 1:parN.nT
                # Write final DOF vector to file
                if i == parN.nT
                    dof = build_dof(state)
                    writedlm("dof-par-$par-trial-$trial-step-$i.csv", dof)
                end
                # Compute network statistics
                Force[i] = network_force(state, state_old, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) # Force vector
                Aggregation_Statistic[i], Bundle_Statistic[i], Parallel_Statistic[i], bundled_directions = bundle_statistic(state, af) # Compute bundle structure statistics
                Bundle_Direction[i] = mean(bundled_directions) # Compute principle bundle direction
                Motor_Position[i] = mean_relative_mp(state)
                Curvature[i] = curvature(af, state, Lxx, Lxy, Lyx, Lyy)
                Motor_Angle[i], TwoFilamentIndex[i] = two_filament_index(mm, state, parA, Lxx, Lxy, Lyx, Lyy)
                # Draw network
                if mod(i-1, parN.plot_interval) == 0
                    draw_network(state, af, xl, mm, parN, parA, Lxx, Lxy, Lyx, Lyy, Force[i], false)
                    savefig("actomyosin-par-$par-trial-$trial-step-$i.pdf") # Save network image
                    dof = build_dof(state)
                    writedlm("dof-par-$par-trial-$trial-step-$i.csv", dof)
                    # draw_bundle_histogram(bundled_directions, Force[i], Lxx, Lyy, par, trial, i) # Draw and save bundle histogram
                end
                # Simulate one time step
                if i != parN.nT # Ensure correct looping sequence
                    af, mm, state = actin_biased_turnover(state, af, mm, parN, parA, Lxx, Lyy) # Implement biased filament turnover
                    af, mm, state = actin_branching_turnover(state, af, mm, parN, parA, Lxx, Lyy) # Implement filament branching
                    af, mm, state = actin_uniform_turnover(state, af, mm, parN, parA, Lxx, Lyy) # Implement uniform filament turnover
                    af, mm, state = actin_treadmilling(state, af, mm, parN, parA) # Implement filament treadmilling
                    af = segment_translations(state, af) # Update periodic BC information
                    xl = intersection_search(state, af, mm) # Update cross-links
                    mm, xl, state = myosin_turnover(state, xl, mm, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) # Myosin binding/unbinding
                end
            end
            # Post-processing
            NetStats = network_statistics(parN, Force, Aggregation_Statistic, Bundle_Statistic, Parallel_Statistic, Bundle_Direction, Motor_Position, Curvature, Motor_Angle, TwoFilamentIndex, Lxx, Lyy)
            draw_statistics(NetStats, parN, parN.nT, par, trial)
        end
    end
    return NetStats
end