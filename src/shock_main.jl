using BenchmarkTools;
using Plots;

""" HW7: Riemann solvers
    For part a, use solver_type = "exact"
    For part b, use solver_type = "approx" """

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
    #solver_type = "approx" # Options: exact/approx
    solver_type = "exact"

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

    if solver_type == "exact"
        riemannsolver = pyimport("riemannsolver") # Credit to Vandenbroucke
        solver = Riemann_solver(γ, riemannsolver)
    elseif solver_type == "approx"
        solver = HLL_solver(γ)
    end

    return benchmark, artificial_viscosity, solver, inputs
end

# ~~~~~~~~~~~~~~~~~ Run  ~~~~~~~~~~~~~~~~~ #
benchmark, artificial_viscosity, solver, inputs = setup()

if benchmark
    # _grid, shock_results, shock_init = @btime run_sim(inputs,
    #                                                     artificial_viscosity)
else
    _grid, shock_results, shock_init = run_sim(inputs,
                                                solver,
                                                artificial_viscosity,
                                                printout=true)
end

shock_anim = to_shock_animation(_grid, shock_results, shock_init, frame_interval=1)
gif(shock_anim, "shock.gif", fps = 15)
