# Generate random variables for thermal motion
# Alex Tam, 09/04/2021

"Generate thermal motion"
function thermal(af, state)
    random = Vector{}();
    for f in af
        nodes = Vector{}(undef, length(state.an[f.index]));
        for i in eachindex(nodes)
            nodes[i] = Vector{}(undef,2); # Pre-allocate
            nodes[i][1] = randn(); nodes[i][2] = randn(); # Store random number
        end
        push!(random, nodes)
    end
    return random
end