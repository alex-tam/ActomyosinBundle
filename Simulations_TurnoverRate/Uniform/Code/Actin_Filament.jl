# Data structures and methods for actin filaments
# Alex Tam, 12/10/2020

"Data structure for actin filament segments"
mutable struct Segment
    index::Int # Integer segment index
    L_eq::Float64 # [μm] Segment equilibrium length
    t::Vector{Vector{Int}} # Translations required by periodic BCs
end

"Mutable data structure for actin filaments"
mutable struct Actin_Filament
    index::Int # Integer filament index
    segments::Vector{Segment} # List of segments on actin filament
    "Inner constructor to generate actin filaments and update State"
    function Actin_Filament(s, ind, parA, Lxx, Lyy, bias::Bool, branch::Bool, nA::Int)
        ###########################################
        ##### New filament at random position #####
        ###########################################
        if (bias == false) && (branch == false) # Generate a filament with random position
            # 1. Generate actin filament segments
            segments = Vector{Segment}(undef, parA.nSeg) # Pre-allocate segments
            for i = 1:parA.nSeg
                segments[i] = Segment(i, parA.LSeg, []) # Create segments
            end
            # 2. Generate random centre position and orientation
            cx = rand() # Dimensionless actin filament centre position (x)
            cy = rand() # Dimensionless actin filament centre position (y)
            angle = 2*pi*rand() # Generate random actin filament orientation
            # 3. Calculate minus end position
            nodes = Vector{}(undef, parA.nSeg+1) # Pre-allocate
            Lf = parA.nSeg*parA.LSeg # Total filament length
            nx = cx - 0.5*cos(angle)*Lf/Lxx # Compute dimensionless x-position of minus end
            ny = cy - 0.5*sin(angle)*Lf/Lyy # Compute dimensionless y-position of minus end
            nodes[1] = Vector{}(undef,2) # Pre-allocate
            nodes[1][1] = nx # Store minus end position (x)
            nodes[1][2] = ny # Store minus end position (y)
            # 4. Calculate other node positions
            for i = 1:parA.nSeg
                nx += segments[i].L_eq*cos(angle)/Lxx # Compute node x-position
                ny += segments[i].L_eq*sin(angle)/Lyy # Compute node y-position
                nodes[i+1] = Vector{}(undef,2) # Pre-allocate
                nodes[i+1][1] = nx # Store node position (x)
                nodes[i+1][2] = ny # Store node position (y)
            end
            push!(s.an, nodes) # Store nodes in state
            # 5. Apply periodic BCs to each segment
            for i = 1:parA.nSeg
                segments[i].t = periodic(s.an[ind][i+1], s.an[ind][i]) # Store required translations
            end
            return new(ind, segments)
        ###########################################
        ##### New filament at biased position #####
        ###########################################
        elseif (bias == true) && (branch == false) # Generate a filament with biased position
            # Specify bias properties
            θ = parA.σ_a*randn() # Filament rotation angle
            dx = parA.σ_p*randn() # Filament translation (x)
            dy = parA.σ_p*randn() # Filament translation (y)
            # 1. Generate actin filament segments
            segments = Vector{Segment}(undef, parA.nSeg) # Pre-allocate segments
            for i = 1:parA.nSeg
                segments[i] = Segment(i, parA.LSeg, []) # Create segments
            end
            # 2. Select reference filament
            ref_index = rand(1:nA) # Random index of an existing filament
            # Find centre point of reference filament
            if mod(parA.nSeg, 2) == 0 # If there is an even number of segments, use central node
                seg::Int = parA.nSeg/2 + 1
                crx = s.an[ref_index][seg][1] # Centre-point of reference filament (x)
                cry = s.an[ref_index][seg][2] # Centre-point of reference filament (y)
            else # If there is an odd number of segments, interpolate
                seg_m::Int = (parA.nSeg+1)/2
                seg_p::Int = (parA.nSeg+1)/2 + 1
                crx = 0.5*(s.an[ref_index][seg_m][1] + s.an[ref_index][seg_p][1]) # Centre-point of reference filament (x)
                cry = 0.5*(s.an[ref_index][seg_m][2] + s.an[ref_index][seg_p][2]) # Centre-point of reference filament (y)
            end 
            # 3. Generate node positions
            nodes = Vector{}(undef, parA.nSeg+1) # Pre-allocate node positions
            # New filament is reference filament rotated about its centre by angle θ
            for i = 1:(parA.nSeg+1)
                nodes[i] = Vector{}(undef,2) # Pre-allocate node position vector
                nodes[i][1] = crx + (s.an[ref_index][i][1] - crx)*cos(θ) - (s.an[ref_index][i][2] - cry)*sin(θ) + dx # Rotate node about (crx, cry) and translate (x)
                nodes[i][2] = cry + (s.an[ref_index][i][1] - crx)*sin(θ) + (s.an[ref_index][i][2] - cry)*cos(θ) + dy # Rotate node about (crx, cry) and translate (y)
            end
            # Reverse polarity if applicable
            rev_test = rand() # Random number [0,1]
            if rev_test <= parA.PRev # If reverse filament polarity
                nodes = reverse(nodes) # Reverse order of nodes
            end
            # Store nodes in state
            push!(s.an, nodes) 
            # 4. Apply periodic BCs to each segment
            for i = 1:parA.nSeg
                segments[i].t = periodic(s.an[ind][i+1], s.an[ind][i]) # Store required translations
            end
            return new(ind, segments)
        #####################################
        ##### New filament by branching #####
        #####################################
        elseif (bias == false) && (branch == true) # Generate a filament by branching
            # 1. Generate actin filament segments
            segments = Vector{Segment}(undef, parA.nSeg) # Pre-allocate segments
            for i = 1:parA.nSeg
                segments[i] = Segment(i, parA.LSeg, []) # Create segments
            end
            # 2. Select reference filament and branchin site
            ind_ref = rand(1:nA) # Random index of an existing filament
            rel_pos = rand() # Branch from random relative position along filaments
            # Locate binding site
            seg = ceil(Int, parA.nSeg*rel_pos) # segment number
            seg_rel_pos = parA.nSeg*rel_pos - seg + 1
            nx = s.an[ind_ref][seg][1] + seg_rel_pos*(s.an[ind_ref][seg+1][1]-s.an[ind_ref][seg][1]) # Compute dimensionless x-position of minus end
            ny = s.an[ind_ref][seg][2] + seg_rel_pos*(s.an[ind_ref][seg+1][2]-s.an[ind_ref][seg][2]) # Compute dimensionless y-position of minus end
            # 3. Generate node positions
            nodes = Vector{}(undef, parA.nSeg+1) # Pre-allocate node positions
            nodes[1] = Vector{}(undef,2) # Pre-allocate minus-end position
            nodes[1][1] = nx # Store minus end position (x)
            nodes[1][2] = ny # Store minus end position (y)
            # Calculate branching angle
            ref_angle = acos((s.an[ind_ref][end][1] - s.an[ind_ref][1][1])/sqrt((s.an[ind_ref][end][1] - s.an[ind_ref][1][1])^2 + (s.an[ind_ref][end][2] - s.an[ind_ref][1][2])^2)) + 0.05*randn();
            left_right = rand() # Random number indicating branching to left or right
            if left_right <= 0.5
                angle = ref_angle + parA.θb/180*pi
            else
                angle = ref_angle - parA.θb/180*pi
            end
            # Calculate other node positions
            for i = 1:parA.nSeg
                nx += segments[i].L_eq*cos(angle)/Lxx # Compute node x-position
                ny += segments[i].L_eq*sin(angle)/Lyy # Compute node y-position
                nodes[i+1] = Vector{}(undef,2) # Pre-allocate
                nodes[i+1][1] = nx # Store node position (x)
                nodes[i+1][2] = ny # Store node position (y)
            end
            push!(s.an, nodes) # Store nodes in state
            # 4. Apply periodic BCs to each segment
            for i = 1:parA.nSeg
                segments[i].t = periodic(s.an[ind][i+1], s.an[ind][i]) # Store required translations
            end
            return new(ind, segments)
        else
            print("Illegal combination for filament turnover.")
            exit()
        end
    end
end

"Obtain dimensional segment lengths"
function get_segment_lengths(f::Actin_Filament, s::State{T}, Lxx, Lxy, Lyx, Lyy) where {T}
    Ls = Vector{T}(undef, length(f.segments)) # Pre-allocate segment lengths
    for seg in f.segments
        mx, my, px, py = get_segment_nodes(f, seg, s, Lxx, Lxy, Lyx, Lyy) # Obtain segment nodes (dimensional, un-translated)
        Ls[seg.index] = sqrt((px - mx)^2 + (py - my)^2) # Compute segment length
    end
    return Ls
end

"Obtain dimensional length of a particular segment"
function get_segment_length(f::Actin_Filament, s::State{T}, seg::Segment, Lxx, Lxy, Lyx, Lyy) where {T}
    mx, my, px, py = get_segment_nodes(f, seg, s, Lxx, Lxy, Lyx, Lyy) # Obtain segment nodes (dimensional, un-translated)
    Ls = sqrt((px - mx)^2 + (py - my)^2) # Compute segment length
    return Ls
end

"Obtain dimensionless, un-translated node positions"
function get_segment_nodes(f::Actin_Filament, seg::Segment, s::State{T}) where {T}
    mx = s.an[f.index][seg.index][1] # Minus end position (x)
    my = s.an[f.index][seg.index][2] # Minus end position (y)
    px = s.an[f.index][seg.index+1][1] # Plus end position (x)
    py = s.an[f.index][seg.index+1][2] # Plus end position (y)
    return mx, my, px, py
end

"Obtain dimensionless, translated node positions"
function get_segment_nodes(f::Actin_Filament, seg::Segment, s::State{T}, t::Vector{Int}) where {T}
    mx = s.an[f.index][seg.index][1] - t[1] # Minus end position (x)
    my = s.an[f.index][seg.index][2] - t[2] # Minus end position (y)
    px = s.an[f.index][seg.index+1][1] - t[1] # Plus end position (x)
    py = s.an[f.index][seg.index+1][2] - t[2] # Plus end position (y)
    return mx, my, px, py
end

"Obtain dimensional, un-translated node positions"
function get_segment_nodes(f::Actin_Filament, seg::Segment, s::State{T}, Lxx, Lxy, Lyx, Lyy) where {T}
    mx = s.an[f.index][seg.index][1]*Lxx + s.an[f.index][seg.index][2]*Lyx # Minus end position (x)
    my = s.an[f.index][seg.index][1]*Lxy + s.an[f.index][seg.index][2]*Lyy # Minus end position (y)
    px = s.an[f.index][seg.index+1][1]*Lxx + s.an[f.index][seg.index+1][2]*Lyx # Plus end position (x)
    py = s.an[f.index][seg.index+1][1]*Lxy + s.an[f.index][seg.index+1][2]*Lyy # Plus end position (y)
    return mx, my, px, py
end

"Obtain dimensional, translated node positions"
function get_segment_nodes(f::Actin_Filament, seg::Segment, s::State{T}, t::Vector{Int}, Lxx, Lxy, Lyx, Lyy) where {T}
    mx = (s.an[f.index][seg.index][1] - t[1])*Lxx + (s.an[f.index][seg.index][2] - t[2])*Lyx # Minus end position (x)
    my = (s.an[f.index][seg.index][1] - t[1])*Lxy + (s.an[f.index][seg.index][2] - t[2])*Lyy # Minus end position (y)
    px = (s.an[f.index][seg.index+1][1] - t[1])*Lxx + (s.an[f.index][seg.index+1][2] - t[2])*Lyx # Plus end position (x)
    py = (s.an[f.index][seg.index+1][1] - t[1])*Lxy + (s.an[f.index][seg.index+1][2] - t[2])*Lyy # Plus end position (y)
    return mx, my, px, py
end