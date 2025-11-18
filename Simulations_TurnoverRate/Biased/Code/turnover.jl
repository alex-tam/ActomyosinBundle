# Implement filament and motor turnover
# Alex Tam, 19/07/2022

"Actin uniform turnover"
function actin_uniform_turnover(s, af, mm, parN, parA, Lxx, Lyy)
    N = length(af) # [-] Number of filaments
    j = 0 # Initialise index counter
    Pt = 1 - exp(-parA.k_uniform*parN.dt) # Turnover probability for any filament
    for i = 1:N # Loop over filaments
        j += 1 # Index (dynamic) of current filament
        random = rand() # Random number [0, 1]
        # Remove filaments
        if random <= Pt
            deleteat!(s.an, j) # Remove DOF from state structure
            deleteat!(af, j) # Remove filament from list
            # Update motors
            for m in mm
                # Mark motors for removal if filament turns over
                if m.f1.index == j
                    m.to_remove = true
                end
                if m.f2.index == j
                    m.to_remove = true
                end
                # Update index of filaments associated with motors
                if m.f1.index > j
                    m.f1.index -= 1 # Reduce index by 1 if necessary
                end
                if m.f2.index > j
                    m.f2.index -= 1 # Reduce index by 1 if necessary
                end
            end
            # Update indices in actin filament list
            for a in eachindex(af)
                af[a].index = a
            end
            j -= 1 # Ensure correct index for next filament
        end
        # Add new filaments
        nNew = parN.nA - length(af) # Number of new filaments (same as number of removed filaments)
        nOld = length(af) # Initial number of filaments for indexing
        for i = 1:nNew
            push!(af, Actin_Filament(s, nOld+i, parA, Lxx, Lyy, false, false, nOld)) # Create new filament and state data
        end
    end
    return af, mm, s
end

"Actin biased turnover"
function actin_biased_turnover(s, af, mm, parN, parA, Lxx, Lyy)
    N = length(af) # [-] Number of filaments
    j = 0 # Initialise index counter
    Pt = 1 - exp(-parA.k_bias*parN.dt) # Turnover probability for any filament
    for i = 1:N # Loop over filaments
        j += 1 # Index (dynamic) of current filament
        random = rand() # Random number [0, 1]
        # Remove filaments
        if random <= Pt
            deleteat!(s.an, j) # Remove DOF from state structure
            deleteat!(af, j) # Remove filament from list
            # Update motors
            for m in mm
                # Mark motors for removal if filament turns over
                if m.f1.index == j
                    m.to_remove = true
                end
                if m.f2.index == j
                    m.to_remove = true
                end
                # Update index of filaments associated with motors
                if m.f1.index > j
                    m.f1.index -= 1 # Reduce index by 1 if necessary
                end
                if m.f2.index > j
                    m.f2.index -= 1 # Reduce index by 1 if necessary
                end
            end
            # Update indices in actin filament list
            for a in eachindex(af)
                af[a].index = a
            end
            j -= 1 # Ensure correct index for next filament
        end
        # Add new filaments
        nNew = parN.nA - length(af) # Number of new filaments (same as number of removed filaments)
        nOld = length(af) # Initial number of filaments for indexing
        for i = 1:nNew
            push!(af, Actin_Filament(s, nOld+i, parA, Lxx, Lyy, true, false, nOld)) # Create new filament and state data
        end
    end
    return af, mm, s
end

"Actin branching turnover"
function actin_branching_turnover(s, af, mm, parN, parA, Lxx, Lyy)
    N = length(af) # [-] Number of filaments
    j = 0 # Initialise index counter
    Pt = 1 - exp(-parA.k_branching*parN.dt) # Turnover probability for any filament
    for i = 1:N # Loop over filaments
        j += 1 # Index (dynamic) of current filament
        random = rand() # Random number [0, 1]
        # Remove filaments
        if random <= Pt
            deleteat!(s.an, j) # Remove DOF from state structure
            deleteat!(af, j) # Remove filament from list
            # Update motors
            for m in mm
                # Mark motors for removal if filament turns over
                if m.f1.index == j
                    m.to_remove = true
                end
                if m.f2.index == j
                    m.to_remove = true
                end
                # Update index of filaments associated with motors
                if m.f1.index > j
                    m.f1.index -= 1 # Reduce index by 1 if necessary
                end
                if m.f2.index > j
                    m.f2.index -= 1 # Reduce index by 1 if necessary
                end
            end
            # Update indices in actin filament list
            for a in eachindex(af)
                af[a].index = a
            end
            j -= 1 # Ensure correct index for next filament
        end
        # Add new filaments
        nNew = parN.nA - length(af) # Number of new filaments (same as number of removed filaments)
        nOld = length(af) # Initial number of filaments for indexing
        for i = 1:nNew
            push!(af, Actin_Filament(s, nOld+i, parA, Lxx, Lyy, false, true, nOld)) # Create new filament and state data
        end
    end
    return af, mm, s
end

"Implement actin filament treadmilling"
function actin_treadmilling(s, af, mm, parN, parA)
    Pt = 1 - exp(-parA.k_t*parN.dt) # Probability of filament treadmilling
    for i = 1:parN.nA # Loop over actin filaments
        random = rand() # Random number [0,1] to simulate treadmilling
        # Apply treadmilling at designated rate
        if random < Pt
            # Add and remove filament segments
            deleteat!(s.an[i], 1) # Remove minus-end segment from filament
            push!(s.an[i], 2*s.an[i][end] - s.an[i][end-1] ) # Add new segment from plus end
            # Update motor relative positions along treadmilled filaments
            for m in mm
                if m.f1.index == i
                    s.mp[m.index][1] -= 1/parA.nSeg # Update motor relative position
                    if s.mp[m.index][1] < 0
                        m.to_remove = true # Remove motor if it was attached to segment now removed
                    end
                elseif m.f2.index == i
                    s.mp[m.index][2] -= 1/parA.nSeg; # Update motor relative position
                    if s.mp[m.index][2] < 0
                        m.to_remove = true # Remove motor if it was attached to segment now removed
                    end
                end
            end
        end
    end
    return af, mm, s
end

"Myosin motor unbinding and unbinding"
function myosin_turnover(s, xl, mm, parN, parA, parM, Lxx, Lxy, Lyx, Lyy)
    N = length(mm) # [-] Number of motors
    j = 0 # Initialise index counter
    for i = 1:N # Loop over motors
        j += 1 # Index (dynamic) of current motor
        random::Float64 = 1 # Pre-allocate for correct scope
        Pt::Float64 = 0 # Pre-allocate for correct scope
        if mm[j].to_remove == false
            random = rand() # Generate random number [0, 1] for current motor
            x1, y1, x2, y2 = get_motor_pos(mm[j], s, parA, Lxx, Lxy, Lyx, Lyy) # Extract dimensional motor positions
            Lm = sqrt((x1-x2)^2 + (y1-y2)^2) # Compute motor length
            off_rate = parM.k_off*exp(parM.k*Lm/parM.F_ref) # Compute off rate of current motor using Bell's Law
            Pt = 1 - exp(-off_rate*parN.dt) # Probability that random turnover will occur at least once in the time step
        end
        # Implement unbinding and immediately replace with new motor
        if any(s.mp[j] .>= 1) || (random <= Pt) || (mm[j].to_remove == true)
            # Remove motor
            deleteat!(s.mp, j) # Remove DOF from state structure
            deleteat!(mm, j) # Remove motor from list
            # Generate replacement motor
            xl_index = rand(1:length(xl)) # Choose random intersection
            push!(mm, Myosin_Motor(s, xl[xl_index], j, Lxx, Lxy, Lyx, Lyy)) # Add motor
            deleteat!(xl, xl_index) # Remove corresponding cross-link
            # Update order of motor indices
            for m = 1:N
                mm[m].index = m
            end
            j -= 1 # Ensure correct index for next motor
        end
    end
    return mm, xl, s
end