function set_shock_ics(grid)
    t = 0
    total_cells = length(grid)
    Q = zeros(Float64, 3, total_cells) # creates an array of zeros with the same shape as x

    # Init velocities
    u = zeros(total_cells)

    # Set boundary condition
    for i in range(1, total_cells, step=1)
        if grid[i] <= 50
            Q[1, i] = 2
        else
            Q[1, i] = 1
        end
        Q[3, i] = 1
    end

    return t, Q
end

function shock_sound_speed(Q, params)
    u = Q[2, :]./Q[1, :]
    ϵₜ = Q[3, :]./Q[1, :]
    ϵₖ = 0.5*(u.^2)
    ϵ = ϵₜ - ϵₖ
    γ = params["γ"]
    return sqrt.(γ*(γ - 1)*ϵ)
end

function shock_pressure_update(Q, dt, dx, istart, iend, params)
    γ = params["γ"]

    # Pressure gradient
    u = Q[2, :]./Q[1, :]
    ϵₜ = Q[3, :]./Q[1, :]
    ϵₖ = 0.5*(u.^2)
    ϵ = ϵₜ - ϵₖ
    p = (γ - 1)*Q[1, :].*ϵ # p = (γ - 1)ρϵ

    # No change in Q[1] (density)
    Q[2, istart:iend] -= (dt/(2*dx))*(p[istart+1:iend+1] - p[istart-1:iend-1])
    Q[3, istart:iend] -= (dt/(2*dx))*(p[istart+1:iend+1].*u[istart+1:iend+1] - p[istart-1:iend-1].*u[istart-1:iend-1])
    return Q
end
