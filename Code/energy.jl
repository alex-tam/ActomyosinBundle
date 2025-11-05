# Actomyosin network energy functional
# Alex Tam, 12/10/2020

"Energy functional"
function energy_functional(x::Vector{T}, s_old::State{Float64}, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) where{T}
    s = build_state(x, af, mm); # Rebuild state from vector input
    energy = zero(T); # Pre-allocate energy
    for f in af
        # energy += energy_actin_thermal(f, s, random, parN, parA, Lxx, Lxy, Lyx, Lyy)
        energy += energy_actin_drag(f, s, s_old, parN, parA, Lxx, Lxy, Lyx, Lyy)
        energy += energy_actin_spring(f, s, parA, Lxx, Lxy, Lyx, Lyy)
        energy += energy_actin_bending(f, s, parA, Lxx, Lxy, Lyx, Lyy)
    end
    for l in xl
        energy += energy_cross_link(l, s, parN, parA, Lxx, Lxy, Lyx, Lyy)
    end
    for m in mm
        energy += energy_myosin_spring(m, s, parA, parM, Lxx, Lxy, Lyx, Lyy)
        energy += energy_myosin_actin(m, s, s_old, parN, parA, parM, Lxx, Lxy, Lyx, Lyy)
    end
    return energy
end

"Energy contribution of thermal forces on actin filaments"
function energy_actin_thermal(f::Actin_Filament, s::State{T}, random, parN, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T) # Pre-allocate energy
    for i in eachindex(s.an[f.index]) # Loop over nodes
        if i == 1 # Minus-end node
            Ls = parA.LSeg
            px = s.an[f.index][i][1]*Lxx + s.an[f.index][i][2]*Lyx 
            py = s.an[f.index][i][1]*Lxy + s.an[f.index][i][2]*Lyy # Dimensional positions
            energy -= sqrt(parA.kb*parA.T*parA.lambda_a*Ls/parN.dt)*( px*random[f.index][i][1] + py*random[f.index][i][2] ) # Thermal energy contribution
        elseif i == length(s.an[f.index]) # Plus-end node
            Ls = parA.LSeg
            px = s.an[f.index][i][1]*Lxx + s.an[f.index][i][2]*Lyx 
            py = s.an[f.index][i][1]*Lxy + s.an[f.index][i][2]*Lyy # Dimensional positions
            energy -= sqrt(parA.kb*parA.T*parA.lambda_a*Ls/parN.dt)*( px*random[f.index][i][1] + py*random[f.index][i][2] ) # Thermal energy contribution
        else # Interior nodees
            Ls1 = parA.LSeg; Ls2 = parA.LSeg
            px = s.an[f.index][i][1]*Lxx + s.an[f.index][i][2]*Lyx 
            py = s.an[f.index][i][1]*Lxy + s.an[f.index][i][2]*Lyy # Dimensional positions
            energy -= sqrt(parA.kb*parA.T*parA.lambda_a*(Ls1+Ls2)/parN.dt)*( px*random[f.index][i][1] + py*random[f.index][i][2] ) # Thermal energy contribution
        end
    end
    return energy
end

"Energy contribution of actin filament spring forces"
function energy_actin_spring(f::Actin_Filament, s::State{T}, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T) # Pre-allocate energy
    for seg in f.segments # Loop over segments
        Ls = get_segment_length(f, s, seg, Lxx, Lxy, Lyx, Lyy) # Physical segment length
        energy += 0.5*parA.k*(Ls - seg.L_eq)^2 # Actin spring energy
    end
    return energy
end

"Energy contribution of drag between actin filaments and the cytoplasm"
function energy_actin_drag(f::Actin_Filament, s::State{T}, s_old::State{Float64}, parN, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T) # Pre-allocate energy
    for i in eachindex(s.an[f.index]) # Loop over nodes
        if i == 1 # Minus-end node
            barl = parA.LSeg/2
            nx = s.an[f.index][i][1] # Dimensionless position (new, x)
            ny = s.an[f.index][i][2] # Dimensionless position (new, y)
            ox = s_old.an[f.index][i][1] # Dimensionless position (old, x)
            oy = s_old.an[f.index][i][2] # Dimensionless position (old, y)
            energy += parA.lambda_a*barl*( ((nx-ox)*Lxx + (ny-oy)*Lyx)^2 + ((nx-ox)*Lxy + (ny-oy)*Lyy)^2 )/(2*parN.dt) # Drag energy contribution
        elseif i == length(s.an[f.index]) # Plus-end node
            barl = parA.LSeg/2
            nx = s.an[f.index][i][1]; ny = s.an[f.index][i][2]; ox = s_old.an[f.index][i][1]; oy = s_old.an[f.index][i][2] # Dimensionless positions
            energy += parA.lambda_a*barl*( ((nx-ox)*Lxx + (ny-oy)*Lyx)^2 + ((nx-ox)*Lxy + (ny-oy)*Lyy)^2 )/(2*parN.dt) # Drag energy contribution
        else # Interior nodees
            barl = parA.LSeg
            nx = s.an[f.index][i][1]; ny = s.an[f.index][i][2]; ox = s_old.an[f.index][i][1]; oy = s_old.an[f.index][i][2] # Dimensionless positions
            energy += parA.lambda_a*barl*( ((nx-ox)*Lxx + (ny-oy)*Lyx)^2 + ((nx-ox)*Lxy + (ny-oy)*Lyy)^2 )/(2*parN.dt) # Drag energy contribution
        end
    end
    return energy
end

"Energy contribution of actin filament bending"
function energy_actin_bending(f::Actin_Filament, s::State{T}, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T) # Pre-allocate energy
    for i = 2:length(s.an[f.index])-1 # Loop over interior nodes
        ddfx, ddfy = second_deriv(f, s, i, parA, Lxx, Lxy, Lyx, Lyy) # Obtain second derivative
        energy += 0.5*parA.kappa*parA.LSeg*(ddfx^2 + ddfy^2) # Add contribution to energy
    end
    return energy
end

"Energy contribution of cross-linker drag"
function energy_cross_link(l::Cross_Link, s::State{T}, parN, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T)
    # Extract cross-linker position (dimensionless, un-translated)
    x1, y1, x2, y2 = get_xl_pos(l, s)
    # Convert to dimensional, translated positions
    lx1 = (x1 - l.t1[1])*Lxx + (y1 - l.t1[2])*Lyx
    ly1 = (x1 - l.t1[1])*Lxy + (y1 - l.t1[2])*Lyy
    lx2 = (x2 - l.t2[1])*Lxx + (y2 - l.t2[2])*Lyx
    ly2 = (x2 - l.t2[1])*Lxy + (y2 - l.t2[2])*Lyy
    # Energy
    energy += parA.lambda_pf*((lx1 - lx2)^2 + (ly1 - ly2)^2)/(2*parN.dt) # Cross-linker drag energy
    return energy
end

"Energy contribution of myosin spring forces"
function energy_myosin_spring(m::Myosin_Motor, s::State{T}, parA, parM, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T)
    # Extract motor binding sites (dimensionless, un-translated)
    x1, y1, x2, y2 = get_motor_pos(m, s, parA, Lxx, Lxy, Lyx, Lyy)
    # Convert to dimensional, translated positions
    mx1 = (x1 - m.t1[1])*Lxx + (y1 - m.t1[2])*Lyx
    my1 = (x1 - m.t1[1])*Lxy + (y1 - m.t1[2])*Lyy
    mx2 = (x2 - m.t2[1])*Lxx + (y2 - m.t2[2])*Lyx
    my2 = (x2 - m.t2[1])*Lxy + (y2 - m.t2[2])*Lyy
    # Energy contribution
    energy += 0.5*parM.k*((mx1 - mx2)^2 + (my1 - my2)^2); # Myosin spring energy
    return energy
end

"Energy contribution of actin--myosin interactions for a single motor"
function energy_myosin_actin(m::Myosin_Motor, s::State{T}, s_old::State{Float64}, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T) # Pre-allocate energy
    Lf = parA.nSeg*parA.LSeg # Filament length
    # Linear component
    energy += parM.Fs/parM.Vm*( Lf*(s.mp[m.index][1] - s_old.mp[m.index][1] ))^2/( 2*parN.dt )
    energy += parM.Fs/parM.Vm*( Lf*(s.mp[m.index][2] - s_old.mp[m.index][2] ))^2/( 2*parN.dt )
    # Constant component
    energy -= parM.Fs*Lf*(s.mp[m.index][1] - s_old.mp[m.index][1])
    energy -= parM.Fs*Lf*(s.mp[m.index][2] - s_old.mp[m.index][2])
    return energy
end

# Obtain dimensional second derivatives z"
function second_deriv(f::Actin_Filament, s::State{T}, i::Int, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    if i == 1 # Minus-end: z'' = 0
        ddfx = 0.0
        ddfy = 0.0
    elseif i == length(s.an[f.index]) # Plus-end: z'' = 0
        ddfx = 0.0 
        ddfy = 0.0
    else # Interior nodes
        # Extract physical positions of plus node
        xp = s.an[f.index][i+1][1]*Lxx + s.an[f.index][i+1][2]*Lyx
        yp = s.an[f.index][i+1][1]*Lxy + s.an[f.index][i+1][2]*Lyy
        # Extract physical positions of centre node
        xc = s.an[f.index][i][1]*Lxx + s.an[f.index][i][2]*Lyx
        yc = s.an[f.index][i][1]*Lxy + s.an[f.index][i][2]*Lyy
        # Extract physical positions of minus node
        xm = s.an[f.index][i-1][1]*Lxx + s.an[f.index][i-1][2]*Lyx
        ym = s.an[f.index][i-1][1]*Lxy + s.an[f.index][i-1][2]*Lyy
        # Compute numerical second derivatives
        ddfx = (xp - 2*xc + xm)/parA.LSeg^2 # Numerical second derivative (x)
        ddfy = (yp - 2*yc + ym)/parA.LSeg^2 # Numerical second derivative (y)
    end
    return ddfx, ddfy
end