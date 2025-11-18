# Data structure and methods for myosin motors
# Alex Tam, 12/10/2020

"Mutable data structure for a myosin motor"
mutable struct Myosin_Motor
    index::Int # Integer motor index
    f1::Actin_Filament # 1st actin filament
    f2::Actin_Filament # 2nd actin filament
    t1::Vector{Int} # 1st filament translation
    t2::Vector{Int} # 2nd filament translation
    to_remove::Bool # Binary variable: true if motor is to be removed (filament turnover)
    "Inner constructor to generate new myosin motors and update State"
    function Myosin_Motor(s::State, l::Cross_Link, ind::Int, Lxx, Lxy, Lyx, Lyy)
        push!(s.mp, get_motor_relative_pos_filament(s, l, Lxx, Lxy, Lyx, Lyy)) # Write relative positions to State
        return new(ind, l.f1, l.f2, l.t1, l.t2, false) # Generate motor
    end
end

"Convert relative position along segments to relative position along filaments"
function get_motor_relative_pos_filament(s::State, l::Cross_Link, Lxx, Lxy, Lyx, Lyy)
    # Compute dimensional segment lengths
    Ls1 = get_segment_lengths(l.f1, s, Lxx, Lxy, Lyx, Lyy) # Filament 1
    Ls2 = get_segment_lengths(l.f2, s, Lxx, Lxy, Lyx, Lyy) # Filament 2
    # Compute relative positions
    mp = Vector{}(undef, 2); # Pre-allocate vector of relative positions
    mp[1] = (sum(Ls1[1:l.seg1_index]) - (1-l.relative_pos1)*Ls1[l.seg1_index])/sum(Ls1)
    mp[2] = (sum(Ls2[1:l.seg2_index]) - (1-l.relative_pos2)*Ls2[l.seg2_index])/sum(Ls2)
    return mp
end

"Convert relative position along filaments to relative position along segments"
function get_motor_relative_pos_segment(m::Myosin_Motor, s::State{T}, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    # 1. Obtain filament lengths (Must be dimensional for non-square domains)
    L1 = parA.nSeg*parA.LSeg
    L2 = parA.nSeg*parA.LSeg
    # 2. Determine dimensional motor positions
    m1 = s.mp[m.index][1]*L1
    m2 = s.mp[m.index][2]*L2
    # 3. Determine index of relevant segments
    L1_cumulative = zero(T) 
    L2_cumulative = zero(T) # Pre-allocate cumulative lengths
    seg1::Int = length(m.f1.segments)
    seg2::Int = length(m.f2.segments) # Pre-allocate
    for i in eachindex(m.f1.segments)
        L1_cumulative += parA.LSeg
        if m1 <= L1_cumulative
            seg1 = i
            break
        end
    end
    for i in eachindex(m.f2.segments)
        L2_cumulative += parA.LSeg
        if m2 <= L2_cumulative
            seg2 = i
            break
        end
    end
    # 4. Determine relative position along segment
    L1_minus = L1_cumulative - parA.LSeg
    L2_minus = L2_cumulative - parA.LSeg
    rel_mp1 = (m1 - L1_minus)/(L1_cumulative - L1_minus)
    rel_mp2 = (m2 - L2_minus)/(L2_cumulative - L2_minus)
    return seg1, seg2, rel_mp1, rel_mp2
end

"Obtain dimensionless, not translated motor positions."
function get_motor_pos(m::Myosin_Motor, s::State{T}, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    # Find relative position along segment
    seg1, seg2, relative_pos1, relative_pos2 = get_motor_relative_pos_segment(m, s, parA, Lxx, Lxy, Lyx, Lyy)
    # Obtain segment nodes position
    m1x, m1y, p1x, p1y = get_segment_nodes(m.f1, m.f1.segments[seg1], s)
    m2x, m2y, p2x, p2y = get_segment_nodes(m.f2, m.f2.segments[seg2], s)
    # Compute motor positions
    x1 = m1x + relative_pos1*(p1x - m1x)
    y1 = m1y + relative_pos1*(p1y - m1y)
    x2 = m2x + relative_pos2*(p2x - m2x)
    y2 = m2y + relative_pos2*(p2y - m2y)
    return x1, y1, x2, y2
end