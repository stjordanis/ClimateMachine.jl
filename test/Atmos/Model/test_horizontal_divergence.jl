#!/usr/bin/env julia --project
using ClimateMachine
ClimateMachine.init(parse_clargs = true)

using ClimateMachine.Atmos
using ClimateMachine.Orientations
using ClimateMachine.ConfigTypes
using ClimateMachine.NumericalFluxes
using ClimateMachine.Diagnostics
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.TurbulenceClosures
using ClimateMachine.SystemSolvers: ManyColumnLU
using ClimateMachine.Mesh.Filters
using ClimateMachine.Mesh.Grids
using ClimateMachine.Mesh.Interpolation
using ClimateMachine.TemperatureProfiles
using ClimateMachine.VariableTemplates
using ClimateMachine.Thermodynamics: air_density, total_energy

using LinearAlgebra
using StaticArrays
using Test

using CLIMAParameters
using CLIMAParameters.Planet: day, planet_radius
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

function init_solid_body_rotation!(problem, bl, state, aux, coords, t)
    FT = eltype(state)

    # initial velocity profile (we need to transform the vector into the Cartesian
    # coordinate system)
    # u_0::FT = 0
    # u_sphere = SVector{3, FT}(u_0, 0, 0)
    # u_init = sphr_to_cart_vec(bl.orientation, u_sphere, aux)
    # e_kin::FT = 0.5 * sum(abs2.(u_init))

    # # Assign state variables
    # state.ρ = aux.ref_state.ρ
    # # state.ρu = u_init
    # state.ρu = u_init .* 0
    # state.ρe = aux.ref_state.ρe + state.ρ * e_kin

    state.ρ = aux.ref_state.ρ
    state.ρu = SVector{3, FT}(0, 0, 0)
    state.ρe = aux.ref_state.ρe

    nothing
end

function config_solid_body_rotation(FT, poly_order, resolution, ref_state)

    # Set up the atmosphere model
    exp_name = "SolidBodyRotation"
    # domain_height::FT = 30e3 # distance between surface and top of atmosphere (m)
    domain_height = FT(50e3)

    model = AtmosModel{FT}(
        AtmosGCMConfigType,
        param_set;
        init_state_prognostic = init_solid_body_rotation!,
        ref_state = ref_state,
        turbulence = ConstantDynamicViscosity(FT(0)),
        moisture = DryModel(),
        source = (Gravity(),),
    )

    config = ClimateMachine.AtmosGCMConfiguration(
        exp_name,
        poly_order,
        resolution,
        domain_height,
        param_set,
        init_solid_body_rotation!;
        model = model,
        numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
    )

    return config
end

function main()
    # Driver configuration parameters
    FT = Float64                             # floating type precision
    poly_order = 3                           # discontinuous Galerkin polynomial order
    n_horz = 12                              # horizontal element number
    n_vert = 6                               # vertical element number
    n_days::FT = 0.05
    timestart::FT = 0                        # start time (s)
    timeend::FT = n_days * day(param_set)    # end time (s)
    @show timeend

    # Set up a reference state for linearization of equations
    temp_profile_ref = DecayingTemperatureProfile{FT}(param_set)
    ref_state = HydrostaticState(temp_profile_ref)

    # Set up driver configuration
    driver_config = config_solid_body_rotation(FT, poly_order, (n_horz, n_vert), ref_state)

    ode_solver_type = ClimateMachine.ExplicitSolverType(
        solver_method = LSRK54CarpenterKennedy,
    )

    CFL = FT(0.01) # target acoustic CFL number

    # time step is computed such that the horizontal acoustic Courant number is CFL
    solver_config = ClimateMachine.SolverConfiguration(
        timestart,
        timeend,
        driver_config,
        Courant_number = CFL,
        ode_solver_type = ode_solver_type,
        CFL_direction = HorizontalDirection(),
        diffdir = HorizontalDirection(),
        fixed_number_of_steps = 2,
    )

    # initialization
    Qinit = similar(solver_config.Q)
    Qinit .= solver_config.Q

    relative_error =
        norm(solver_config.Q .- Qinit) /
        norm(Qinit)
    @info "Relative error = $relative_error"

    # Set up diagnostics
    dgn_config = config_diagnostics(FT, driver_config)

    cb_print_step = GenericCallbacks.EveryXSimulationSteps(1) do
        @show getsteps(solver_config.solver)
        relative_error =
        norm(solver_config.Q .- Qinit) /
        norm(Qinit)
        @info "Relative error = $relative_error"
        t = gettime(solver_config.solver)
        @info "t = $t"
        nothing
    end

    # Run the model
    result = ClimateMachine.invoke!(
        solver_config;
        diagnostics_config = dgn_config,
        user_callbacks = (cb_print_step,),
        check_euclidean_distance = false,
    )
    t_end = gettime(solver_config.solver)
    @info "t_end = $t_end"

    relative_error =
        norm(solver_config.Q .- Qinit) /
        norm(Qinit)
    @info "Relative error = $relative_error"
    return solver_config, Qinit
end

function config_diagnostics(FT, driver_config)
    interval = "0.5shours" # chosen to allow diagnostics every 30 simulated minutes

    _planet_radius = FT(planet_radius(param_set))

    info = driver_config.config_info
    boundaries = [
        FT(-90.0) FT(-180.0) _planet_radius
        FT(90.0) FT(180.0) FT(_planet_radius + info.domain_height)
    ]
    resolution = (FT(2), FT(2), FT(1000)) # in (deg, deg, m)
    interpol = ClimateMachine.InterpolationConfiguration(
        driver_config,
        boundaries,
        resolution,
    )

    dgngrp = setup_atmos_default_diagnostics(
        AtmosGCMConfigType(),
        interval,
        driver_config.name,
        interpol = interpol,
    )

    return ClimateMachine.DiagnosticsConfiguration([dgngrp])
end

solver_config, Qinit = main();
nothing
