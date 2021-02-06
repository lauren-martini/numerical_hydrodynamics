using BenchmarkTools;
using Plots;

include("run_sim.jl")
include("utils.jl")
include("flow.jl")
include("shock.jl")
include("plot_anim.jl")

# ~~~~~~~~~~~ Simulation settings ~~~~~~~~~~~ #
benchmark = false
#case = "flow"
case = "shock"

# ~~~~~~~~~~~ Simulation inputs ~~~~~~~~~~~ #
domain = [0, 100]
N = 500 # num grid points/cells
ghosts = 2 # num ghost cells

if case == "flow"
    end_time = 100 # length of sim
    sound_speed(Q, params) = 1 # Speed of sound - constant for this isothermal case
elseif case == "shock"
    end_time = 40
    sound_speed = shock_sound_speed
end

# Courant number - see CFL condition (u*dt/dx <= C)
C = 0.5
# Adiabatic exponent - constant
γ = 1.4

inputs = Dict("domain"=>domain,
              "N"=>N,
              "ghosts"=>ghosts,
              "end_time"=>end_time,
              "C"=>C,
              "sound_speed"=>sound_speed,
              "γ"=>γ)

# ~~~~~~~~~~~~~~~~~ Run  ~~~~~~~~~~~~~~~~~ #
if case == "flow"
    if benchmark
        _grid, flow_results, flow_init = @btime run_sim(inputs,
                                                set_flow_ics,
                                                flow_pressure_update)
    else
        _grid, flow_results, flow_init = run_sim(inputs,
                                               set_flow_ics,
                                               flow_pressure_update, printout=true)
    end

    flow_anim = to_flow_animation(_grid, flow_results, flow_init, frame_interval=10)
    gif(flow_anim, "flow.gif", fps = 15)

elseif case == "shock"
    if benchmark
        _grid, shock_results, shock_init = @btime run_sim(inputs,
                                               set_shock_ics,
                                               shock_pressure_update)
    else
        _grid, shock_results, shock_init = run_sim(inputs,
                                               set_shock_ics,
                                               shock_pressure_update, printout=true)
    end

    shock_anim = to_shock_animation(_grid, shock_results, shock_init, frame_interval=10)
    gif(shock_anim, "shock.gif", fps = 15)
end

# plot(_grid, flow_results[1, :], ylims=(0.95, 1.35), color="grey", label="Q")
# #plot!(size=(900, 400))
# plot(_grid, flow_init[1, :], label="Initial Condition", line=[:dash], color="blue")
# xlabel!("x")
# ylabel!("Q(x, t)")

# plot(_grid, shock_results[1, :], ylims=(0.95, 2.1), color="grey", label="Q")
# #plot!(size=(900, 400))
# plot!(_grid, shock_init[1, :], label="Initial Condition", line=[:dash], color="blue")
# xlabel!("x")
# ylabel!("Q(x, t)")
