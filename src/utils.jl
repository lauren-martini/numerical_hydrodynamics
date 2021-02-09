function calc_soundspeed(γ, p, ρ)
    return  sqrt.(γ*p./ρ)
end

function Riemann_solver(γ, riemannsolver)
    """ Factory for the riemann solver so that
    its output matches that of the HLL_solver. """
    solver = riemannsolver.RiemannSolver(γ).solve
    function riemann_solver(ρₗ, uₗ, pₗ, ρᵣ, uᵣ, pᵣ)
        ρ_sol, u_sol, p_sol, _tmp = solver(ρₗ, uₗ, pₗ,
                                           ρᵣ, uᵣ, pᵣ)
        flux = [0.0, 0.0, 0.0]
        flux[1] = ρ_sol*u_sol
        flux[2] = (ρ_sol*u_sol)^2 + p_sol
        flux[3] = (p_sol*γ ./ (γ - 1) .+ 0.5*ρ_sol*u_sol.^2)*u_sol
        return flux
    end
end

function HLL_solver(γ)
    """ Factory for hll_solvers so that their
    input matches that in Reimannsolver.py"""

    function hll_solver(ρₗ, uₗ, pₗ, ρᵣ, uᵣ, pᵣ)
        # Compute sound speeds at interfaces
        cₛₗ = calc_soundspeed(γ, pₗ, ρₗ)
        cₛᵣ = calc_soundspeed(γ, pᵣ, ρᵣ)

        # Compute signal speeds
        sₗ = min(uₗ - cₛₗ, uᵣ - cₛᵣ)
        sᵣ = max(uₗ + cₛₗ, uᵣ + cₛᵣ)

        # Compute flux
        if sₗ <= 0 <= sᵣ
            Fₗ = [ρₗ*uₗ,
                 (ρₗ*uₗ)^2 + pₗ,
                 (pₗ*γ / (γ - 1) + 0.5*ρₗ*uₗ^2)*uₗ]
            Fᵣ = [ρᵣ*uᵣ,
                  (ρᵣ*uᵣ)^2 + pᵣ,
                  (pᵣ*γ / (γ - 1) + 0.5*ρᵣ*uᵣ^2)*uᵣ]
            F_hll = [0.0, 0.0, 0.0]
            F_hll[1] = (sᵣ*Fₗ[1] - sₗ*Fᵣ[1] .+ sₗ*sᵣ*(ρᵣ - ρₗ))/(sᵣ - sₗ)
            F_hll[2] = (sᵣ*Fₗ[2] - sₗ*Fᵣ[2] .+ sₗ*sᵣ*(uᵣ - uₗ))/(sᵣ - sₗ)
            F_hll[3] = (sᵣ*Fₗ[3] - sₗ*Fᵣ[3] .+ sₗ*sᵣ*(pᵣ - pₗ))/(sᵣ - sₗ)
        elseif 0 <= sₗ # Fl
            F_hll = [ρₗ*uₗ,
                    (ρₗ*uₗ)^2 + pₗ,
                    (pₗ*γ / (γ - 1) + 0.5*ρₗ*uₗ^2)*uₗ]
        elseif 0 >= sᵣ # Fr
            F_hll = [ρᵣ*uᵣ,
                    (ρᵣ*uᵣ)^2 + pᵣ,
                    (pᵣ*γ / (γ - 1) + 0.5*ρᵣ*uᵣ^2)*uᵣ]
        else
            print("Something went wrong. Check HLL solver.")
        end

        return F_hll
    end
end

function gen_grid(x, domain, ghosts, dx, total_cells)
    grid = zero(x)
    grid[1] = domain[1] - (dx/2)*(ghosts - 1)
    for i in range(2, total_cells, step=1)
        grid[i] = grid[i - 1] + dx
    end
    return grid
end
