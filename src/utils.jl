# function flux_step(Q, istart, iend)
#     flux = zero(Q)
#     for i in range(istart, iend, step=1)
#         # Velocity across the boundary
#         u_boundary = 0.5((Q[2, i]/Q[1, i]) + (Q[2, i-1]/Q[1, i-1]))
#         # u_boundary = calc_u_boundary(Q, i)
#
#         if u_boundary >= 0
#             flux[:, i] = Q[:, i-1] * u_boundary
#         else
#             flux[:, i] = Q[:, i] * u_boundary
#         end
#     end
#     return flux
# end

function calc_soundspeed(γ, p, ρ)
    return  sqrt.(γ*p./ρ)
end

function Riemann_solver(solver)
    #todo
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
                 (pₗ*uₗ)^2 + pₗ,
                 (pₗ*γ / (γ - 1) + 0.5*ρₗ*uₗ^2)*uₗ]
            Fᵣ = [ρᵣ*uᵣ,
                  (pᵣ*uᵣ)^2 + pᵣ,
                  (pᵣ*γ / (γ - 1) + 0.5*ρᵣ*uᵣ^2)*uᵣ]
            F_hll = (sᵣ*Fₗ - sₗ*Fᵣ .+ sₗ*sᵣ*(uᵣ - uₗ))/(sᵣ - sₗ)
            ρ_sol = F_hll[1]
            u_sol = F_hll[2]
            p_sol = F_hll[3]
        elseif 0 <= sₗ # Fl
            ρ_sol = ρₗ*uₗ
            u_sol = (pₗ*uₗ)^2 + pₗ
            p_sol = (pₗ*γ / (γ - 1) + 0.5*ρₗ*uₗ^2)*uₗ
        elseif 0 >= sᵣ # Fr
            ρ_sol = ρᵣ*uᵣ
            u_sol = (pᵣ*uᵣ)^2 + pᵣ
            p_sol = (pᵣ*γ / (γ - 1) + 0.5*ρᵣ*uᵣ^2)*uᵣ
        else
            print("Something went wrong. Check HLL solver.")
        end

        return ρ_sol, u_sol, p_sol, nothing
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
