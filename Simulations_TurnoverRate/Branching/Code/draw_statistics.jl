# Produce plots of actomyosin network statistics
# Alex Tam, 25/07/2022

"Main function to draw network statistics"
function draw_statistics(NetStats::Network_Statistics, parN, time_step, par, trial)
    draw_stress(NetStats.Stress, parN, time_step, par, trial)
    draw_mean_normal_stress(NetStats.Mean_Normal_Stress, parN, time_step, par, trial)
end

"Draw time-dependent stress components"
function draw_stress(Stress, parN, time_step, par, trial)
    plot() # Clear existing plots
    times = (0:time_step-1).*parN.dt # [s] Vector of times
    σxx = Vector{Float64}() # [pN/um] Pre-allocate stress component vector (σxx)
    σxy = Vector{Float64}() # [pN/um] Pre-allocate stress component vector (σxy)
    σyx = Vector{Float64}() # [pN/um] Pre-allocate stress component vector (σyx)
    σyy = Vector{Float64}() # [pN/um] Pre-allocate stress component vector (σyy)
    for i = 1:time_step # Loop over time steps
        push!(σxx, Stress[i][1])
        push!(σxy, Stress[i][2])
        push!(σyx, Stress[i][3])
        push!(σyy, Stress[i][4])
    end
    # Plot stress components
    plot(times, σxx, color = "orange", linewidth = 1, alpha = 0.25, label=L"\sigma_{xx}", xlabel = L"$t$", ylabel = L"\mathrm{Stress}\;\mathrm{(pN/}\mu\mathrm{m)}", legend=:bottomleft, margin = 5mm)
    plot!(times, σxy, color = "purple", linewidth = 1, alpha = 0.25, label=L"\sigma_{xy}")
    plot!(times, σyx, color = "green", linewidth = 1, alpha = 0.25, label=L"\sigma_{yx}")
    plot!(times, σyy, color = "blue", linewidth = 1, alpha = 0.25, label=L"\sigma_{yy}")
    # Plot stress components (moving averages)
    plot!(times, movingaverage(σxx, 10), color = "orange", linewidth = 1, linestyle = :dash, label = false)
    plot!(times, movingaverage(σxy, 10), color = "purple", linewidth = 1, linestyle = :dash, label = false)
    plot!(times, movingaverage(σyx, 10), color = "green", linewidth = 1, linestyle = :dash, label = false)
    plot!(times, movingaverage(σyy, 10), color = "blue", linewidth = 1, linestyle = :dash, label = false)
    savefig("stress-par-$par-trial-$trial.pdf") # Save image
end

"Draw time-dependent mean normal stress"
function draw_mean_normal_stress(Mean_Normal_Stress, parN, time_step, par, trial)
    plot() # Clear existing plots
    times = (0:time_step-1).*parN.dt # [s] Vector of times
    plot(times, Mean_Normal_Stress, color=:black, linewidth = 1, alpha = 0.25, label=false, xlabel = L"$t$", ylabel = L"\mathrm{Stress}\;\mathrm{(pN/}\mu\mathrm{m)}", legend=:bottomleft, margin = 5mm)
    plot!(times, movingaverage(Mean_Normal_Stress, 10), color=:black, linewidth = 1, linestyle=:dash, label=false)
    savefig("mean_normal_stress-par-$par-trial-$trial.pdf") # Save image
end

"Draw histogram of bundle directions"
function draw_bundle_histogram(bundled_directions, force, Lxx, Lyy, par, trial, i)
    stress = [ force[1]/Lyy force[2]/Lyy ; force[3]/Lxx force[4]/Lxx ] # Compute stress matrix
    evals = eigvals(stress) # Eigenvalues of stress matrix
    evecs = eigvecs(stress) # Eigenvectors of stress matrix
    θ1 = atan(evecs[1,2], evecs[1,1]) # Principal angle corresponding to first eigenvalue, [-π,π]
    θ2 = atan(evecs[2,2], evecs[2,1]) # Principal angle corresponding to second eigenvalue, [-π,π]
    hist = fit(Histogram, bundled_directions, -π:π/20:π)
    histogram(bundled_directions, xlabel = L"$\theta$", bins = range(-π, π, step = π/20), leg = false; init = 0)
    # Plot principal directions (blue = contraction, orange = expansion)
    max_pstress = max(abs(evals[1]), abs(evals[2])) # Maximum principal stress for normalisation
    max_hist = maximum(hist.weights) # Maximum amplitude of histogram
    if evals[1] <= 0
        plot!([θ1, θ1], [0, max_hist*abs(evals[1])/max_pstress], color = "blue", linewidth = 5)
    else
        plot!([θ1, θ1], [0, max_hist*abs(evals[1])/max_pstress], color = "orange", linewidth = 5)
    end
    if evals[2] <= 0
        plot!([θ2, θ2], [0, max_hist*abs(evals[2])/max_pstress], color = "blue", linewidth = 5)
    else
        plot!([θ2, θ2], [0, max_hist*abs(evals[2])/max_pstress], color = "orange", linewidth = 5)
    end
    savefig("bundle_hist-par-$par-trial-$trial-step-$i.pdf") # Save bundle direction histogram
    writedlm("bundle_directions-par-$par-trial-$trial-step-$i.csv", bundled_directions) # Save bundle direction data to file
    writedlm("stress_matrix-par-$par-trial-$trial-step-$i.csv", stress) # Save stress matrix to file
end
