# Initial conditions
function set_flow_ics(grid)
    t = 0
    total_cells = length(grid)
    Q = zeros(Float64, 2, total_cells) # creates an array of zeros with the same shape as x

    # Init velocities
    u = zeros(total_cells)

    # Set boundary condition
    for i in range(1, total_cells, step=1)
        Q[1, i] = 1 + 0.3*exp(-((grid[i] - 50)^2)/10)
        Q[2, i] = Q[1, i]*u[i]
    end

    return t, Q
end

function flow_pressure_update(Q, dt, dx, istart, iend, params)
    # Pressure gradient
    cₛ = params["cₛ"]

    # No change in Q[1] (density)
    dp = (cₛ^2) * (Q[1, istart+1:iend+1] - Q[1, istart-1:iend-1])./(2*dx)
    Q[2, istart:iend] -= dt*dp

    return Q
end
