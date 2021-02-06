# Initial conditions
function set_ics(grid)
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


function run_sim(inputs; printout=false)
    # Unpack simulation inputs
    domain = inputs["domain"]
    N = inputs["N"]
    ghosts = inputs["ghosts"]
    end_time = inputs["end_time"]
    cₛ = inputs["cₛ"]

    total_cells = N + ghosts
    x = collect(range(domain[1], domain[2], length=total_cells))
    dx = (domain[2] - domain[1])/N
    istart = Int(ghosts/2) + 1
    iend = N + Int(ghosts/2)

    # Generate grid
    grid = gen_grid(x, domain, ghosts, dx, total_cells)

    t, Q = set_ics(grid) # case-specific initial conditions
    Q_init = copy(Q)
    results = zeros(Float64, total_cells, 1) # storage for Q at each timepoint

    timestep = 1
    while t < end_time

        # calculate time_step
        dt = C*minimum(dx./(cₛ .+ abs.(Q[2,istart:iend]./Q[1,istart:iend])))

        # ~~~ ADVECTION ~~~ #

        # calculate fluxes
        fluxes = flux_step(Q, istart, iend)

        # Apply advection
        #Q[1, istart:iend] -= (dt/dx)*(fluxes[1, istart+1:iend+1] - fluxes[1, istart:iend])
        Q[:, istart:iend] -= (dt/dx)*(fluxes[:, istart+1:iend+1] - fluxes[:, istart:iend])

        # ~~~ SOURCE TERMS ~~~ #

        # Update Pressure
        # No change in Q[1] (density)
        dp = (cₛ^2) * (Q[1, istart+1:iend+1] - Q[1, istart-1:iend-1])./(2*dx)
        Q[2, istart:iend] -= dt*dp

        # Boundary conditions
        Q[1, 1] = Q[1, 2] # rho_0 = rho_1
        Q[2, 1] = -Q[2, 2]
        Q[1, end] = Q[1, end-1] # rho_N+1 = rho_N
        Q[2, end] = -Q[2, end-1]

        if timestep % 10 == 0
            results = [results Q[1, :]]
        end
        t += dt
        timestep += 1

        if printout
            print("\r t = $t    dt = $(dt)")
        end
    end

    return grid, results, Q_init
end


function to_flow_animation(grid, res, init; frame_interval=10)

    num_res = size(res)[2]
    println()
    anim = @animate for i in 1:num_res
        plot(grid, res[:, i], ylims=(0.95, 1.35), color="grey", label="Density", legend=:topleft)
        #plot!(size=(900, 400))
        plot!(grid, init[1, :], label="Initial Condition", line=[:dash], color="blue")
        xlabel!("x")
        ylabel!("Q(x, t)")
        title!("Frame $i out of $num_res")

        #println()
        print("\rAnimating frame $i out of $num_res")
    end every frame_interval

    return anim
end
