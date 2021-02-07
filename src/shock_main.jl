using BenchmarkTools;
using Plots;

try
    using PyCall;
catch
    println("Adding PyCall. This can take a few minutes.\n")
    Pkg.add("PyCall")
    using PyCall;
end
scriptdir = @__DIR__
pushfirst!(PyVector(pyimport("sys")."path"), scriptdir)

include("shock.jl")
include("utils.jl")

function setup()
    # ~~~~~~~~~~~ Simulation settings ~~~~~~~~~~~ #
    benchmark = false
    artificial_viscosity = false
    riemannsolver = pyimport("riemannsolver") # Credit to Vandenbroucke

    # ~~~~~~~~~~~ Simulation inputs ~~~~~~~~~~~ #
    domain = [-.5, .5]
    N = 500 # num grid points/cells
    ghosts = 2 # num ghost cells
    end_time = 0.2

    # Courant number - see CFL condition (u*dt/dx <= C)
    C = 0.3
    # Adiabatic exponent - constant
    γ = 1.4

    inputs = Dict("domain"=>domain,
                  "N"=>N,
                  "ghosts"=>ghosts,
                  "end_time"=>end_time,
                  "C"=>C,
                  "γ"=>γ)
    return benchmark, artificial_viscosity, riemannsolver, inputs
end

# ~~~~~~~~~~~~~~~~~ Run  ~~~~~~~~~~~~~~~~~ #
benchmark, artificial_viscosity, riemannsolver, inputs = setup()

if benchmark
    # _grid, shock_results, shock_init = @btime run_sim(inputs,
    #                                                     artificial_viscosity)
else
    _grid, shock_results, shock_init = run_sim(inputs,
                                                riemannsolver,
                                                artificial_viscosity,
                                                printout=true)
end

shock_anim = to_shock_animation(_grid, shock_results, shock_init, frame_interval=1)
gif(shock_anim, "shock.gif", fps = 15)
