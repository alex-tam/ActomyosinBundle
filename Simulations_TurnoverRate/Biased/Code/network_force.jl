# Compute network force components
# Alex Tam, 12/10/2020

"Compute the contractile/expansive forces in the network"
function network_force(s, s_old, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy)
    x::Vector{Real} = build_dof(s) # Build DOF vector
    # Compute the network force (+ve expansive, -ve contractile)
    fxx = -ForwardDiff.derivative(a -> energy_functional(x, s_old, af, xl, mm, random, parN, parA, parM, a, Lxy, Lyx, Lyy), Lxx::Real)
    fxy = -ForwardDiff.derivative(a -> energy_functional(x, s_old, af, xl, mm, random, parN, parA, parM, Lxx, a, Lyx, Lyy), Lxy::Real)
    fyx = -ForwardDiff.derivative(a -> energy_functional(x, s_old, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, a, Lyy), Lyx::Real)
    fyy = -ForwardDiff.derivative(a -> energy_functional(x, s_old, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, a), Lyy::Real)
    return [fxx, fxy, fyx, fyy]
end