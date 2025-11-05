# Apply dimensionless periodic boundary conditions
# Alex Tam, 12/10/2020

"Apply dimensionless periodic boundary conditions to segments and motors"
function periodic(pos1, pos2)
    translations = Vector{Vector{Int}}(undef, 5) # Pre-allocate vector of translation vectors
    translations[1] = [0,0] # Ensure that the identity translation occurs first
    translations[2] = [floor(pos1[1]), floor(pos1[2])] # Store 1st possible combination
    translations[3] = [floor(pos1[1]), floor(pos2[2])] # Store 2nd possible combination
    translations[4] = [floor(pos2[1]), floor(pos1[2])] # Store 3rd possible combination
    translations[5] = [floor(pos2[1]), floor(pos2[2])] # Store 4th possible combination
    translations = unique(translations) # Delete repeated translations
    return translations
end

"Update all filament segment translations"
function segment_translations(s, af)
    for f in af
        for seg in f.segments
            seg.t = periodic(s.an[f.index][seg.index+1], s.an[f.index][seg.index]) # Apply periodic BCs
        end
    end
    return af
end