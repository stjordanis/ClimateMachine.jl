# # [Implicit-Explicit (IMEX) Additive Runge-Kutta Timestepping](@id Single-rate-IMEXARK-Timestepping)

# Lore ipsum

include(joinpath(
    @__DIR__,
    "../../../../../",
    "tutorials/Numerics/TimeStepping/tutorial_acousticwave_config.jl",
))

FT = Float64
timeend = FT(3600)

ode_solver =
    ClimateMachine.IMEXSolverType(solver_method = ARK2ImplicitExplicitMidpoint)
C = FT(0.2)
cfl_direction = HorizontalDirection()
run_acousticwave(ode_solver, C, cfl_direction, timeend)
