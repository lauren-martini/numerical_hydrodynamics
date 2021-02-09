function set_ics(grid)
    t = 0
    total_cells = length(grid)
    Q = zeros(Float64, 3, total_cells) # creates an array of zeros with the same shape as x

    # Init velocities
    u = zeros(total_cells)

    # Set boundary condition
    for i in range(1, total_cells, step=1)
        if grid[i] <= 0
            Q[1, i] = 2
        else
            Q[1, i] = 1
        end
        ϵₜ = 1
        Q[2, i] = Q[1, i]*u[i]
        Q[3, i] = Q[1, i]*ϵₜ
    end

    return t, Q
end


function calc_av(Q, u, ζ, istart, iend)
    Π = zero(u)
    for i in range(istart, iend, step=1)
        if u[i + 1] <= u[i - 1]
            Π[i] = 0.25*(ζ^2)*((u[i+1] - u[i-1]).^2).*Q[1, i]
        end
    end
    return Π
end


function run_sim(inputs, solver, artificial_viscosity; printout=false)
    # Unpack simulation inputs
    domain = inputs["domain"]
    N = inputs["N"]
    ghosts = inputs["ghosts"]
    end_time = inputs["end_time"]
    C = inputs["C"]
    γ = inputs["γ"]

    if artificial_viscosity
        ζ = 3.0
    end

    total_cells = N + ghosts
    x = collect(range(domain[1], domain[2], length=total_cells))
    dx = (domain[2] - domain[1])/N
    istart = Int(ghosts/2) + 1
    iend = N + Int(ghosts/2)

    # Generate grid
    grid = gen_grid(x, domain, ghosts, dx, total_cells)

    t, Q = set_ics(grid) # case-specific initial conditions
    Q_init = copy(Q)
    density = zeros(Float64, total_cells, 1) # storage for Q1 at each timepoint
    velocity = zeros(Float64, total_cells, 1)
    energy = zeros(Float64, total_cells, 1)
    pressure = zeros(Float64, total_cells, 1)

    timestep = 1
    while t < end_time

        # Update Pressure
        ρ = Q[1, :]     # Density
        u = Q[2, :]./ρ  # Velocity
        ϵₜ = Q[3, :]./ρ  # Total energy
        ϵₖ = 0.5*(u.^2)  # Kinetic energy
        ϵ = ϵₜ - ϵₖ       # Internal energy

        p = (γ - 1)*ρ.*ϵ # Pressure = (γ - 1)ρϵ
        if artificial_viscosity
            Π = calc_av(Q, u, ζ, istart, iend)
            p += Π
        end

        # calculate time_step
        if any(x->x<0, ρ) || any(x->x<0, p)
            println("\nwarning, negative sound speed")
        end
        cₛ = calc_soundspeed(γ, p, ρ)
        dt = C*minimum(dx./(cₛ[istart:iend] .+ abs.(u[istart:iend])))

        # ~~~ ADVECTION ~~~ #

        # calculate fluxes
        fluxes = zero(Q)
        for i in range(istart, iend+1, step=1)
            fluxes[1, i], fluxes[2, i], fluxes[3, i] = solver(ρ[i-1], u[i-1],
                                                              p[i-1], ρ[i],
                                                              u[i], p[i])
        end

        # Apply advection
        Q[:, istart:iend] -= (dt/dx)*(fluxes[:, istart+1:iend+1] - fluxes[:, istart:iend])

        # Boundary conditions
        Q[1, 1] = Q[1, 2] # rho_0 = rho_1
        Q[2, 1] = -Q[2, 2]
        Q[3, 1] = Q[3, 2]
        Q[1, end] = Q[1, end-1] # rho_N+1 = rho_N
        Q[2, end] = -Q[2, end-1]
        Q[3, end] = Q[3, end-1]

        # Save data
        if timestep % 10 == 0
            density = [density Q[1, :]]
            velocity = [velocity Q[2, :]./Q[1, :]]
            energy = [energy Q[3, :]./Q[1, :]]
            pressure = [pressure p]
        end

        t += dt
        timestep += 1

        if printout
            print("\r t = $t    dt = $(dt)")
            #println("t = $t      dt = $(dt)")
        end
    end

    return grid, [density, velocity, energy, pressure], Q_init
end


function to_shock_animation(grid, res, init; frame_interval=10)

    num_res = size(res[1])[2]
    println("\nPrepping plot for animation.\n")
    anim = @animate for i in 1:num_res
        #plot(grid, res[:, i], ylims=(0.95, 2.1), color="grey", label="Density", legend=:topright)
        #plot!(size=(900, 400))

        # Density
        density_plot = plot(grid, res[1][:, i], label="ρ")
        plot!(grid, init[1, :], line=[:dash], color="blue", label="IC")
        #xlabel!("x")
        ylabel!("Density")
        title!("Frame $i out of $num_res")

        # Velocity
        velocity_plot = plot(grid, res[2][:, i], legend=false)
        ylabel!("Velocity")

        energy_plot = plot(grid, res[3][:, i], legend=false)
        ylabel!("Energy")

        pressure_plot = plot(grid, res[4][:, i], legend=false)
        ylabel!("Pressure")

        plot(density_plot,
                velocity_plot,
                energy_plot,
                pressure_plot, layout=(4, 1))

        print("\rAnimating frame $i out of $num_res")
    end every frame_interval

    return anim
end
