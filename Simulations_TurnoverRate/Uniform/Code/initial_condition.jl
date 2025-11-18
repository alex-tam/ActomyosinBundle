# Generate initial conditions
# Alex Tam, 12/10/2020

"Initial condition for actin filaments"
function actin_ic(s::State, parN, parA, Lxx, Lyy)
    af = Vector{Actin_Filament}(undef, parN.nA)
    for i in eachindex(af)
        af[i] = Actin_Filament(s, i, parA, Lxx, Lyy, false, false, parN.nA) # Generate initial actin filament
    end
    return af, s
end

"Initial condition for myosin motors"
function myosin_ic(s::State, mm, parN, xl, Lxx, Lxy, Lyx, Lyy)
    xl_indices = sample(1:length(xl), parN.nM, replace = false) # Sample intersections without replacement
    for i = 1:parN.nM
        push!(mm, Myosin_Motor(s, xl[xl_indices[i]], i, Lxx, Lxy, Lyx, Lyy)) # Generate myosin motor
    end
    deleteat!(xl, sort(xl_indices)) # Remove cross-links at intersections with a motor
    return mm, xl, s
end