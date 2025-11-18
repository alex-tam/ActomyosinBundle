# Data structure and method for cross-links
# Alex Tam, 12/10/2020

"Data structure for a cross-linker"
struct Cross_Link
    f1::Actin_Filament # 1st cross-linked actin filament
    f2::Actin_Filament # 2nd cross-linked actin filament
    seg1_index::Int # Index of relevant segment on 1st filament
    seg2_index::Int # Index of relevant segment on 2nd filament
    relative_pos1::Float64 # Relative position along the 1st segment
    relative_pos2::Float64 # Relative position along the 2nd segment
    t1::Vector{Int} # Translation of the 1st segment
    t2::Vector{Int} # Translation of the 2nd segment
end

"Obtain dimensionless cross-linker position"
function get_xl_pos(l::Cross_Link, s::State)
    # Find actin segment nodes
    m1x, m1y, p1x, p1y = get_segment_nodes(l.f1, l.f1.segments[l.seg1_index], s)
    m2x, m2y, p2x, p2y = get_segment_nodes(l.f2, l.f2.segments[l.seg2_index], s)
    # Compute the cross-link site
    x1 = m1x + l.relative_pos1*(p1x - m1x)
    y1 = m1y + l.relative_pos1*(p1y - m1y)
    x2 = m2x + l.relative_pos2*(p2x - m2x)
    y2 = m2y + l.relative_pos2*(p2y - m2y)
    return x1, y1, x2, y2
end