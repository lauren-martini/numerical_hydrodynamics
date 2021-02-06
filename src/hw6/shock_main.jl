using BenchmarkTools;
using Plots;

include("shock.jl")
include("utils.jl")

# ~~~~~~~~~~~ Simulation settings ~~~~~~~~~~~ #
benchmark = true
artificial_viscosity = false

"""Yes, there is a noticable decrease in oscillations when
artificial viscosity is turned on. """

# ~~~~~~~~~~~ Simulation inputs ~~~~~~~~~~~ #
domain = [0, 100]
N = 500 # num grid points/cells
ghosts = 2 # num ghost cells
end_time = 40

# Courant number - see CFL condition (u*dt/dx <= C)
C = 0.5
# Adiabatic exponent - constant
γ = 1.4

inputs = Dict("domain"=>domain,
              "N"=>N,
              "ghosts"=>ghosts,
              "end_time"=>end_time,
              "C"=>C,
              "γ"=>γ)

# ~~~~~~~~~~~~~~~~~ Run  ~~~~~~~~~~~~~~~~~ #
if benchmark
    _grid, shock_results, shock_init = @btime run_sim(inputs,
                                                        artificial_viscosity)
else
    _grid, shock_results, shock_init = run_sim(inputs,
                                                artificial_viscosity,
                                                printout=true)
end

shock_anim = to_shock_animation(_grid, shock_results, shock_init, frame_interval=1)
gif(shock_anim, "shock.gif", fps = 15)
