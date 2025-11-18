# Minimise energy functional to compute solutions
# Alex Tam, 12/10/2020

"Minimise energy functional to solve model"
function optimise_network(s_old, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy)
    dof = build_dof(s_old) # Construct DOF vector for initial guess
    @time result = optimize(x -> energy_functional(x, s_old, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy), (G, x) -> energy_gradient!(G, x, s_old, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy), dof, LBFGS(m = 20, alphaguess = LineSearches.InitialHagerZhang(α0 = 1.0)), Optim.Options(x_abstol = parN.xTol, f_abstol = parN.fTol, g_abstol = parN.gTol) )
    if Optim.converged(result) == 0 # Check convergence
        @printf("Warning: Optimisation did not converge.\n")
    end
    return Optim.minimizer(result) # Return DOF vector that minimises energy
end