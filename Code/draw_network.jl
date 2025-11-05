# Produce actomyosin network plots
# Alex Tam, 19/07/2022

"Draw the network configuration"
function draw_network(s, af, xl, mm, parN, parA, Lxx, Lxy, Lyx, Lyy, force, p_stress::Bool)
    plot() # Clear existing plots
    # Actin Filaments
    fx = zeros(Float64, parA.nSeg+1, length(af))
    fy = zeros(Float64, parA.nSeg+1, length(af)) # Pre-allocate un-translated filament node positions
    n1x = Vector{Float64}()
    n1y = Vector{Float64}()
    n2x = Vector{Float64}()
    n2y = Vector{Float64}() # Pre-allocate translated node positions (non-plus end segments)
    p1x = Vector{Float64}()
    p1y = Vector{Float64}() 
    p2x = Vector{Float64}()
    p2y = Vector{Float64}() # Pre-allocate translated node positions (plus end segments)
    for f in af
        nx = Vector{}()
        ny = Vector{}() # Pre-allocate vectors of current filament nodes
        for seg in f.segments
            for j in eachindex(seg.t)
                if j == 1
                    mx, my, px, py = get_segment_nodes(f, seg, s, seg.t[j], Lxx, Lxy, Lyx, Lyy) # Get segment nodes
                    push!(nx, mx) 
                    push!(nx, px) # Store unprojected nodes of a segment (x)
                    push!(ny, my) 
                    push!(ny, py) # Store unprojected nodes of a segment (y)
                else
                    mx, my, px, py = get_segment_nodes(f, seg, s, seg.t[j], Lxx, Lxy, Lyx, Lyy) # Get segment nodes
                    if seg.index == length(f.segments)
                        push!(p1x, mx)
                        push!(p2x, px)
                        push!(p1y, my)
                        push!(p2y, py) # Store nodes (plus-end)
                    else
                        push!(n1x, mx)
                        push!(n2x, px)
                        push!(n1y, my)
                        push!(n2y, py) # Store nodes (non-plus-end)
                    end
                end
            end
        end
        nx = unique(nx) # Remove duplicates
        ny = unique(ny) # Remove duplicates
        fx[:, f.index] = nx # Store filament nodes in matrix
        fy[:, f.index] = ny # Store filament nodes in matrix
    end
    plot!(fx, fy, legend = false, color = "red", aspect_ratio=:equal, xlabel = L"$x$", ylabel = L"$y$", xlims = [0, parN.lxx], ylims = [0, parN.lyy]) # Plot un-translated segments
    scatter!(fx[1:end-1,:], fy[1:end-1,:], color = "red") # Plot nodes as points
    scatter!(fx[end,:], fy[end,:], color = "black") # Plot plus end nodes
    plot!(transpose(hcat(n1x, n2x)), transpose(hcat(n1y, n2y)), color = "red") # Plot translated segments (non-plus-end)
    plot!(transpose(hcat(p1x, p2x)), transpose(hcat(p1y, p2y)), color = "red") # Plot translated segments (plus end)
    scatter!([n1x, n2x], [n1y, n2y], color = "red") # Plot translated nodes (non-plus-ends)
    scatter!([p1x], [p1y], color = "red") # Plot translated nodes
    scatter!([p2x], [p2y], color = "black") # Plot translated nodes (plus ends)
    # Cross-Links
    # Pre-allocate outside loop to plot all  cross-links at once
    l1x = Vector{Float64}()
    l1y = Vector{Float64}() 
    l2x = Vector{Float64}()
    l2y = Vector{Float64}()
    for l in xl
        # Obtain cross-link positions (dimensionless, un-translated)
        x1, y1, x2, y2 = get_xl_pos(l, s)
        # Convert to dimensional, translated positions
        lx1 = (x1 - l.t1[1])*Lxx + (y1 - l.t1[2])*Lyx
        ly1 = (x1 - l.t1[1])*Lxy + (y1 - l.t1[2])*Lyy
        lx2 = (x2 - l.t2[1])*Lxx + (y2 - l.t2[2])*Lyx
        ly2 = (x2 - l.t2[1])*Lxy + (y2 - l.t2[2])*Lyy
        push!(l1x, lx1)
        push!(l1y, ly1)
        push!(l2x, lx2)
        push!(l2y, ly2)
    end
    # scatter!([l1x, l2x], [l1y, l2y], color = "pink") # Plot cross-link sites as points
    # plot!(transpose(hcat(l1x, l2x)), transpose(hcat(l1y, l2y)), color = "pink") # Draw cross-linkers
    # Myosin motors
    # Pre-allocate outside loop to plot all motors at once
    m1x = Vector{Float64}()
    tm1x = Vector{Float64}()
    m1y = Vector{Float64}()
    tm1y = Vector{Float64}()
    m2x = Vector{Float64}()
    tm2x = Vector{Float64}()
    m2y = Vector{Float64}()
    tm2y = Vector{Float64}()
    for m in mm
        x1, y1, x2, y2 = get_motor_pos(m, s, parA, Lxx, Lxy, Lyx, Lyy); # Obtain motor positions (dimensionless, un-translated)
        # Convert to dimensional, translated positions
        mx1 = (x1 - m.t1[1])*Lxx + (y1 - m.t1[2])*Lyx
        my1 = (x1 - m.t1[1])*Lxy + (y1 - m.t1[2])*Lyy
        mx2 = (x2 - m.t2[1])*Lxx + (y2 - m.t2[2])*Lyx
        my2 = (x2 - m.t2[1])*Lxy + (y2 - m.t2[2])*Lyy
        # Compute translations for visualisation only
        myosin_translations = periodic([mx1/parN.lxx, my1/parN.lyy], [mx2/parN.lxx, my2/parN.lyy]);
        for j in eachindex(myosin_translations)
            if j == 1
                push!(m1x, mx1)
                push!(m1y, my1)
                push!(m2x, mx2)
                push!(m2y, my2)
            else
                t = myosin_translations[j]
                # Store dimensional, translated positions
                push!(tm1x, mx1 - (t[1]*Lxx + t[2]*Lyx))
                push!(tm1y, my1 - (t[1]*Lxy + t[2]*Lyy))
                push!(tm2x, mx2 - (t[1]*Lxx + t[2]*Lyx))
                push!(tm2y, my2 - (t[1]*Lxy + t[2]*Lyy))
            end
        end
    end
    scatter!([m1x, m2x], [m1y, m2y], color = "blue") # Plot binding sites as points
    scatter!([tm1x, tm2x], [tm1y, tm2y], color = "blue") # Plot binding sites as points
    plot!(transpose(hcat(m1x, m2x)), transpose(hcat(m1y, m2y)), color = "blue") # Draw motors
    plot!(transpose(hcat(tm1x, tm2x)), transpose(hcat(tm1y, tm2y)), color = "blue") # Draw projected motors
    # Optional: Draw principal stress vectors on network plot
    if p_stress == true
        stress = [ force[1]/Lyy force[2]/Lyy ; force[3]/Lxx force[4]/Lxx ] # Compute stress matrix
        evals = eigvals(stress) # Eigenvalues of stress matrix
        evecs = eigvecs(stress) # Eigenvectors of stress matrix
        # Add principal stress vectors to plot
        if evals[1] >= 0
            plot!([parN.lxx/2, parN.lxx/2 + evals[1]*evecs[1,1]], [parN.lyy/2 , parN.lyy/2 + evals[1]*evecs[1,2]], arrow = :arrow, color = "orange")
        else
            plot!([parN.lxx/2, parN.lxx/2 + abs(evals[1])*evecs[1,1]], [parN.lyy/2 , parN.lyy/2 + abs(evals[1])*evecs[1,2]], arrow = :arrow, color = "blue")
        end
        if evals[2] >= 0
            plot!([parN.lxx/2, parN.lxx/2 + evals[2]*evecs[2,1]], [parN.lyy/2 , parN.lyy/2 + evals[2]*evecs[2,2]], arrow = :arrow, color = "orange")
        else
            plot!([parN.lxx/2, parN.lxx/2 + abs(evals[2])*evecs[2,1]], [parN.lyy/2 , parN.lyy/2 + abs(evals[2])*evecs[2,2]], arrow = :arrow, color = "blue")
        end
    end
end