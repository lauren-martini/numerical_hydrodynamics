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

function gen_grid(x, domain, ghosts, dx, total_cells)
    grid = zero(x)
    grid[1] = domain[1] - (dx/2)*(ghosts - 1)
    for i in range(2, total_cells, step=1)
        grid[i] = grid[i - 1] + dx
    end
    return grid
end
