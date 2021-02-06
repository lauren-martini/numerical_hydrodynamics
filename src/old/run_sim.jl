function run_sim(inputs, case_ics, update_pressure; printout=false)
    # Unpack simulation inputs
    domain = inputs["domain"]
    N = inputs["N"]
    ghosts = inputs["ghosts"]
    end_time = inputs["end_time"]
    sound_speed = inputs["sound_speed"]

    total_cells = N + ghosts
    x = collect(range(domain[1], domain[2], length=total_cells))
    dx = (domain[2] - domain[1])/N
    istart = Int(ghosts/2) + 1
    iend = N + Int(ghosts/2)

    # Generate grid
    grid = zero(x)
    grid[1] = domain[1] - (dx/2)*(ghosts - 1)
    for i in range(2, total_cells, step=1)
        grid[i] = grid[i - 1] + dx
    end

    t, Q = case_ics(grid) # case-specific initial conditions
    Q_init = copy(Q)
    results = zeros(Float64, total_cells, 1) # storage for Q at each timepoint

    while t < end_time

        # calculate time_step
        cₛ = sound_speed(Q, inputs)
        dt = C*minimum(dx./(cₛ[istart:iend] .+ abs.(Q[2,istart:iend]./Q[1,istart:iend])))

        # ~~~ ADVECTION ~~~ #

        # calculate fluxes
        fluxes = flux_step(Q, istart, iend)

        # Apply advection
        #Q[1, istart:iend] -= (dt/dx)*(fluxes[1, istart+1:iend+1] - fluxes[1, istart:iend])
        Q[:, istart:iend] -= (dt/dx)*(fluxes[:, istart+1:iend+1] - fluxes[:, istart:iend])

        # ~~~ SOURCE TERMS ~~~ #

        # Update Pressure
        Q = update_pressure(Q, dt, dx, istart, iend, inputs)

        # Boundary conditions
        # Q[1, 1] = Q[1, 2] # rho_0 = rho_1
        # Q[2, 1] = -Q[2, 2]
        # Q[1, end] = Q[1, end-1] # rho_N+1 = rho_N
        # Q[2, end] = -Q[2, end-1]

        results = [results Q[1, :]]
        t += dt

        if printout
            print("\r t = $t    dt = $(dt)")
        end
    end

    return grid, results, Q_init
end
