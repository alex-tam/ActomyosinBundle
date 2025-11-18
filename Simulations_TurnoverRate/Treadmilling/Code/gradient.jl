# Analytical gradient for the energy functional
# Alex Tam, 01/07/2021

function energy_gradient!(G::Vector{Float64}, x::Vector{Float64}, s_old::State{Float64}, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy)
    s = build_state(x, af, mm); # Rebuild state from vector input
    grad_vec = zeros(Float64, size(x)); grad_state = build_state(grad_vec, af, mm); # Place gradient terms in a state structure
    for f in af
        # gradient_actin_thermal!(grad_state, f, s, random, parN, parA, Lxx, Lxy, Lyx, Lyy)
        gradient_actin_drag!(grad_state, f, s, s_old, parN, parA, Lxx, Lxy, Lyx, Lyy)
        gradient_actin_spring!(grad_state, f, s, parA, Lxx, Lxy, Lyx, Lyy)
        gradient_actin_bending!(grad_state, f, s, parA, Lxx, Lxy, Lyx, Lyy)
    end
    for l in xl
        gradient_cross_link!(grad_state, l, s, parN, parA, Lxx, Lxy, Lyx, Lyy)
    end
    for m in mm
        gradient_myosin_spring!(grad_state, m, s, parA, parM, Lxx, Lxy, Lyx, Lyy)
        gradient_myosin_actin!(grad_state, m, s, s_old, parN, parA, parM, Lxx, Lxy, Lyx, Lyy)
    end
    grad_vec += build_dof(grad_state)
    for j in eachindex(grad_vec)
        G[j] = grad_vec[j]
    end
    return G
end

"Gradient contribution for thermal motion"
function gradient_actin_thermal!(grad_state, f::Actin_Filament, s::State{T}, random, parN, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    for i in eachindex(s.an[f.index]) # Loop over nodes
        if i == 1 # Minus-end node
            barl = parA.LSeg/2
            px = s.an[f.index][i][1]*Lxx; py = s.an[f.index][i][2]*Lyy # Dimensional positions
            grad_state.an[f.index][i][1] += -sqrt(2*parA.kb*parA.T*parA.lambda_a*barl/parN.dt)*random[f.index][i][1]*Lxx # Thermal forces contribution
            grad_state.an[f.index][i][2] += -sqrt(2*parA.kb*parA.T*parA.lambda_a*barl/parN.dt)*random[f.index][i][2]*Lyy # Thermal forces contribution
        elseif i == length(s.an[f.index]) # Plus-end node
            px = s.an[f.index][i][1]*Lxx; py = s.an[f.index][i][2]*Lyy # Dimensional positions
            barl = parA.LSeg/2
            grad_state.an[f.index][i][1] += -sqrt(2*parA.kb*parA.T*parA.lambda_a*barl/parN.dt)*random[f.index][i][1]*Lxx # Thermal forces contribution
            grad_state.an[f.index][i][2] += -sqrt(2*parA.kb*parA.T*parA.lambda_a*barl/parN.dt)*random[f.index][i][2]*Lyy # Thermal forces contribution
        else # Interior nodes
            px = s.an[f.index][i][1]*Lxx; py = s.an[f.index][i][2]*Lyy # Dimensional positions
            barl = parA.LSeg
            grad_state.an[f.index][i][1] += -sqrt(2*parA.kb*parA.T*parA.lambda_a*barl/parN.dt)*random[f.index][i][1]*Lxx # Thermal forces contribution
            grad_state.an[f.index][i][2] += -sqrt(2*parA.kb*parA.T*parA.lambda_a*barl/parN.dt)*random[f.index][i][2]*Lyy # Thermal forces contribution
        end
    end
    return build_dof(grad_state)
end

"Gradient contribution for drag between actin filaments and the cytoplasm"
function gradient_actin_drag!(grad_state, f::Actin_Filament, s::State{T}, s_old::State{Float64}, parN, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    for i in eachindex(s.an[f.index]) # Loop over nodes
        if i == 1 # Minus-end node
            barl = parA.LSeg/2
            xn = s.an[f.index][i][1] # Segment minus node current position (x)
            yn = s.an[f.index][i][2] # Segment minus node current position (y)
            xo = s_old.an[f.index][i][1] # Segment minus node old position (x)
            yo = s_old.an[f.index][i][2] # Segment minus node old position (y)
            grad_state.an[f.index][i][1] += parA.lambda_a*barl*(xn-xo)*Lxx^2/(parN.dt) # Add drag for minus node (x)
            grad_state.an[f.index][i][2] += parA.lambda_a*barl*(yn-yo)*Lyy^2/(parN.dt) # Add drag for minus node (y)
        elseif i == length(s.an[f.index]) # Plus-end node
            barl = parA.LSeg/2
            xn = s.an[f.index][i][1] # Segment minus node current position (x)
            yn = s.an[f.index][i][2] # Segment minus node current position (y)
            xo = s_old.an[f.index][i][1] # Segment minus node old position (x)
            yo = s_old.an[f.index][i][2] # Segment minus node old position (y)
            grad_state.an[f.index][i][1] += parA.lambda_a*barl*(xn-xo)*Lxx^2/(parN.dt) # Add drag for minus node (x)
            grad_state.an[f.index][i][2] += parA.lambda_a*barl*(yn-yo)*Lyy^2/(parN.dt) # Add drag for minus node (y)
        else # Interior nodes
            barl = parA.LSeg
            xn = s.an[f.index][i][1] # Segment minus node current position (x)
            yn = s.an[f.index][i][2] # Segment minus node current position (y)
            xo = s_old.an[f.index][i][1] # Segment minus node old position (x)
            yo = s_old.an[f.index][i][2] # Segment minus node old position (y)
            grad_state.an[f.index][i][1] += parA.lambda_a*barl*(xn-xo)*Lxx^2/(parN.dt) # Add drag for minus node (x)
            grad_state.an[f.index][i][2] += parA.lambda_a*barl*(yn-yo)*Lyy^2/(parN.dt) # Add drag for minus node (y)
        end
    end
    return grad_state
end

"Gradient contribution of actin filament spring forces"
function gradient_actin_spring!(grad_state, f::Actin_Filament, s::State{T}, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    for i in eachindex(s.an[f.index]) # Loop over nodes
        if i == 1 # Minus-end node
            Ls = get_segment_length(f, s, f.segments[1], Lxx, Lxy, Lyx, Lyy) # Physical segment length
            grad_state.an[f.index][i][1] += -Lxx^2*(s.an[f.index][2][1] - s.an[f.index][1][1])*parA.k*(Ls - f.segments[1].L_eq)/Ls
            grad_state.an[f.index][i][2] += -Lyy^2*(s.an[f.index][2][2] - s.an[f.index][1][2])*parA.k*(Ls - f.segments[1].L_eq)/Ls 
        elseif i == length(s.an[f.index]) # Plus-end node
            Ls = get_segment_length(f, s, f.segments[end], Lxx, Lxy, Lyx, Lyy) # Physical segment length
            grad_state.an[f.index][i][1] += Lxx^2*(s.an[f.index][end][1] - s.an[f.index][end-1][1])*parA.k*(Ls - f.segments[end].L_eq)/Ls 
            grad_state.an[f.index][i][2] += Lyy^2*(s.an[f.index][end][2] - s.an[f.index][end-1][2])*parA.k*(Ls - f.segments[end].L_eq)/Ls
        else # Interior nodes
            Ls1 = get_segment_length(f, s, f.segments[i-1], Lxx, Lxy, Lyx, Lyy) # Physical segment length
            Ls2 = get_segment_length(f, s, f.segments[i], Lxx, Lxy, Lyx, Lyy) # Physical segment length
            grad_state.an[f.index][i][1] += -Lxx^2*(s.an[f.index][i+1][1] - s.an[f.index][i][1])*parA.k*(Ls2 - f.segments[i].L_eq)/Ls2 + Lxx^2*(s.an[f.index][i][1] - s.an[f.index][i-1][1])*parA.k*(Ls1 - f.segments[i-1].L_eq)/Ls1 # (x)
            grad_state.an[f.index][i][2] += -Lyy^2*(s.an[f.index][i+1][2] - s.an[f.index][i][2])*parA.k*(Ls2 - f.segments[i].L_eq)/Ls2 + Lyy^2*(s.an[f.index][i][2] - s.an[f.index][i-1][2])*parA.k*(Ls1 - f.segments[i-1].L_eq)/Ls1 # (y)
        end
    end
    return grad_state
end

"Gradient contribution of actin filament bending"
function gradient_actin_bending!(grad_state,  f::Actin_Filament, s::State{T}, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    for i in eachindex(s.an[f.index]) # Loop over interior nodes
        if i == 1 # Minus-end node
            ddspx, ddspy = second_deriv(f, s, i+1, parA, Lxx, Lxy, Lyx, Lyy)
            grad_state.an[f.index][i][1] += parA.kappa*Lxx/parA.LSeg*ddspx
            grad_state.an[f.index][i][2] += parA.kappa*Lyy/parA.LSeg*ddspy
        elseif i == length(s.an[f.index]) # Plus-end node
            ddsmx, ddsmy = second_deriv(f, s, i-1, parA, Lxx, Lxy, Lyx, Lyy)
            grad_state.an[f.index][i][1] += parA.kappa*Lxx/parA.LSeg*ddsmx
            grad_state.an[f.index][i][2] += parA.kappa*Lyy/parA.LSeg*ddsmy
        else # Interior nodes
            # Obtain second derivatives
            ddsmx, ddsmy = second_deriv(f, s, i-1, parA, Lxx, Lxy, Lyx, Lyy)
            ddscx, ddscy = second_deriv(f, s, i, parA, Lxx, Lxy, Lyx, Lyy)
            ddspx, ddspy = second_deriv(f, s, i+1, parA, Lxx, Lxy, Lyx, Lyy)
            grad_state.an[f.index][i][1] += parA.kappa*Lxx/parA.LSeg*(ddspx - 2*ddscx + ddsmx) # x-component
            grad_state.an[f.index][i][2] += parA.kappa*Lyy/parA.LSeg*(ddspy - 2*ddscy + ddsmy) # y-component
        end
    end
    return grad_state
end

"Gradient contribution of protein friction"
function gradient_cross_link!(grad_state, l::Cross_Link, s::State{T}, parN, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    x1, y1, x2, y2 = get_xl_pos(l::Cross_Link, s::State) # Dimensionless
    # First filament
    grad_state.an[l.f1.index][l.seg1_index][1] += parA.lambda_pf*Lxx^2/(parN.dt)*(x1-l.t1[1]-x2+l.t2[1])*(1-l.relative_pos1)
    grad_state.an[l.f1.index][l.seg1_index][2] += parA.lambda_pf*Lyy^2/(parN.dt)*(y1-l.t1[2]-y2+l.t2[2])*(1-l.relative_pos1)
    grad_state.an[l.f1.index][l.seg1_index+1][1] += parA.lambda_pf*Lxx^2/(parN.dt)*(x1-l.t1[1]-x2+l.t2[1])*(l.relative_pos1)
    grad_state.an[l.f1.index][l.seg1_index+1][2] += parA.lambda_pf*Lyy^2/(parN.dt)*(y1-l.t1[2]-y2+l.t2[2])*(l.relative_pos1)
    # Second filament
    grad_state.an[l.f2.index][l.seg2_index][1] += -parA.lambda_pf*Lxx^2/(parN.dt)*(x1-l.t1[1]-x2+l.t2[1])*(1-l.relative_pos2)
    grad_state.an[l.f2.index][l.seg2_index][2] += -parA.lambda_pf*Lyy^2/(parN.dt)*(y1-l.t1[2]-y2+l.t2[2])*(1-l.relative_pos2)
    grad_state.an[l.f2.index][l.seg2_index+1][1] += -parA.lambda_pf*Lxx^2/(parN.dt)*(x1-l.t1[1]-x2+l.t2[1])*(l.relative_pos2)
    grad_state.an[l.f2.index][l.seg2_index+1][2] += -parA.lambda_pf*Lyy^2/(parN.dt)*(y1-l.t1[2]-y2+l.t2[2])*(l.relative_pos2)
    return grad_state
end

"Gradient contribution of myosin spring forces"
function gradient_myosin_spring!(grad_state, m, s::State{T}, parA, parM, Lxx, Lxy, Lyx, Lyy) where {T}
    # Effect on filament positions
    s1, s2, m1, m2 = get_motor_relative_pos_segment(m::Myosin_Motor, s::State{T}, parA, Lxx, Lxy, Lyx, Lyy)
    x1, y1, x2, y2 = get_motor_pos(m, s, parA, Lxx, Lxy, Lyx, Lyy) # Dimensionless, not translated
    # First filament
    grad_state.an[m.f1.index][s1][1] += parM.k*Lxx^2*(x1-m.t1[1]-x2+m.t2[1])*(1-m1)
    grad_state.an[m.f1.index][s1][2] += parM.k*Lyy^2*(y1-m.t1[2]-y2+m.t2[2])*(1-m1)
    grad_state.an[m.f1.index][s1+1][1] += parM.k*Lxx^2*(x1-m.t1[1]-x2+m.t2[1])*m1
    grad_state.an[m.f1.index][s1+1][2] += parM.k*Lyy^2*(y1-m.t1[2]-y2+m.t2[2])*m1
    # Second filament
    grad_state.an[m.f2.index][s2][1] += -parM.k*Lxx^2*(x1-m.t1[1]-x2+m.t2[1])*(1-m2)
    grad_state.an[m.f2.index][s2][2] += -parM.k*Lyy^2*(y1-m.t1[2]-y2+m.t2[2])*(1-m2)
    grad_state.an[m.f2.index][s2+1][1] += -parM.k*Lxx^2*(x1-m.t1[1]-x2+m.t2[1])*m2
    grad_state.an[m.f2.index][s2+1][2] += -parM.k*Lyy^2*(y1-m.t1[2]-y2+m.t2[2])*m2
    # Motor position
    grad_state.mp[m.index][1] += parM.k*(x1-m.t1[1]-x2+m.t2[1])*Lxx^2*(s.an[m.f1.index][s1+1][1] - s.an[m.f1.index][s1][1])*parA.nSeg + parM.k*(y1-m.t1[2]-y2+m.t2[2])*Lyy^2*(s.an[m.f1.index][s1+1][2] - s.an[m.f1.index][s1][2])*parA.nSeg
    grad_state.mp[m.index][2] += -parM.k*(x1-m.t1[1]-x2+m.t2[1])*Lxx^2*(s.an[m.f2.index][s2+1][1] - s.an[m.f2.index][s2][1])*parA.nSeg - parM.k*(y1-m.t1[2]-y2+m.t2[2])*Lyy^2*(s.an[m.f2.index][s2+1][2] - s.an[m.f2.index][s2][2])*parA.nSeg
    return grad_state
end

"Gradient contribution of actin--myosin interactions for a single motor"
function gradient_myosin_actin!(grad_state, m::Myosin_Motor, s::State{T}, s_old::State{Float64}, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) where {T}
    Lf = parA.nSeg*parA.LSeg # Filament length
    # Gradient
    grad_state.mp[m.index][1] += parM.Fs/parM.Vm*Lf*(s.mp[m.index][1] - s_old.mp[m.index][1])/parN.dt - parM.Fs*Lf
    grad_state.mp[m.index][2] += parM.Fs/parM.Vm*Lf*(s.mp[m.index][2] - s_old.mp[m.index][2])/parN.dt - parM.Fs*Lf
    return grad_state
end