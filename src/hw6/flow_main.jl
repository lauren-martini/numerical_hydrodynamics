using BenchmarkTools;
using Plots;

include("flow.jl")
include("utils.jl")

# ~~~~~~~~~~~ Simulation settings ~~~~~~~~~~~ #
benchmark = false

# ~~~~~~~~~~~ Simulation inputs ~~~~~~~~~~~ #
domain = [0, 100]
N = 500 # num grid points/cells
ghosts = 2 # num ghost cells
end_time = 100 # length of sim

# Sound speed - constant for isothermal case
cₛ = 1
# Courant number - see CFL condition (u*dt/dx <= C)
C = 0.5

inputs = Dict("domain"=>domain,
              "N"=>N,
              "ghosts"=>ghosts,
              "end_time"=>end_time,
              "C"=>C,
              "cₛ"=>cₛ)

# ~~~~~~~~~~~~~~~~~ Run  ~~~~~~~~~~~~~~~~~ #
if benchmark
    _grid, flow_results, flow_init = @btime run_sim(inputs)
else
    _grid, flow_results, flow_init = run_sim(inputs, printout=true)
end

flow_anim = to_flow_animation(_grid, flow_results, flow_init, frame_interval=1)
gif(flow_anim, "flow.gif", fps = 15)
