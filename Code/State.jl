# Data structure and methods for degrees of freedom
# Alex Tam, 12/10/2020

"Mutable data structure of degrees of freedom"
mutable struct State{T}
    an::Vector{Vector{Vector{T}}} # Dimensionless positions (x, y) of actin filament nodes
    mp::Vector{Vector{T}} # Relative positions (m1, m2) of myosin motors along attached filaments
end

"Convert the state structure to a vector for use in the optimiser"
function build_dof(s::State{Float64})
    # Type specification is required here
    an = Vector{Float64}() # Pre-allocate filament positions
    mp = Vector{Float64}() # Pre-allocate motor positions
    # Actin filament node positions
    for i = 1:length(s.an) # Loop over filaments
        for j = 1:length(s.an[i]) # Loop over nodes
            push!(an, s.an[i][j][1]) # Store node x-position
            push!(an, s.an[i][j][2]) # Store node y-position
        end
    end
    # Myosin motor relative positions
    for i = 1:length(s.mp) # Loop over motors
        push!(mp, s.mp[i][1]) # Store motor position on 1st filament
        push!(mp, s.mp[i][2]) # Store motor position on 2nd filament
    end
    return vcat(an, mp)
end

"Convert the state vector to a data structure for easier access"
function build_state(x::Vector{T}, af, mm) where {T}
    an = Vector{Vector{Vector{T}}}(undef, length(af)) # Pre-allocate filament positions
    mp = Vector{Vector{T}}(undef, length(mm)) # Pre-allocate motor positions
    # Actin filament node positions
    dof_count = 0 # Initialise DOF counter
    for i in eachindex(af) # Loop over filaments
        nodes = Vector{Vector{T}}(undef, length(af[i].segments)+1) # Pre-allocate all nodes on an current filament
        for j = 1:(length(af[i].segments)+1) # Loop over nodes
            nodes[j] = Vector{T}(undef,2) # Pre-allocate current node position
            nodes[j][1] = x[dof_count+1] # Obtain node positions
            nodes[j][2] = x[dof_count+2] # Obtain node positions
            dof_count += 2; # Increment counter
        end
        an[i] = nodes; # Construct vector of vectors of vectors
    end
    # Myosin motor relative positions
    for i in eachindex(mm) # Loop over motors
        mp[i] = Vector{T}(undef, 2) # Pre-allocate current motor position
        mp[i][1] = x[dof_count + (i-1)*2 + 1] # Obtain motor positions
        mp[i][2] = x[dof_count + (i-1)*2 + 2] # Obtain motor positions
    end
    return State{T}(an, mp)
end