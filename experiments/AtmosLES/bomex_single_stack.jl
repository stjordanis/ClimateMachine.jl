include("bomex_model.jl")

function main()
    # add a command line argument to specify the kind of surface flux
    # TODO: this will move to the future namelist functionality
    bomex_args = ArgParseSettings(autofix_names = true)
    add_arg_group!(bomex_args, "BOMEX")
    @add_arg_table! bomex_args begin
        "--moisture-model"
        help = "specify cloud condensate model"
        metavar = "equilibrium|nonequilibrium"
        arg_type = String
        default = "equilibrium"
        "--surface-flux"
        help = "specify surface flux for energy and moisture"
        metavar = "prescribed|bulk"
        arg_type = String
        default = "prescribed"
    end

    cl_args =
        ClimateMachine.init(parse_clargs = true, custom_clargs = bomex_args)

    surface_flux = cl_args["surface_flux"]
    moisture_model = cl_args["moisture_model"]

    FT = Float64
    config_type = SingleStackConfigType

    # DG polynomial order
    N = 1

    # Prescribe domain parameters
    nelem_vert = 50
    zmax = FT(3000)

    t0 = FT(0)

    # For a full-run, please set the timeend to 3600*6 seconds
    # For the test we set this to == 30 minutes
    timeend = FT(1800)
    #timeend = FT(3600 * 6)
    CFLmax = FT(0.90)

    # Choose default IMEX solver
    ode_solver_type = ClimateMachine.IMEXSolverType()

    model = bomex_model(
        FT,
        config_type,
        zmax,
        surface_flux;
        moisture_model = moisture_model,
    )
    ics = model.problem.init_state_prognostic

    # Assemble configuration
    driver_config = ClimateMachine.SingleStackConfiguration(
        "BOMEX_SINGLE_STACK",
        N,
        nelem_vert,
        zmax,
        param_set,
        model;
        solver_type = ode_solver_type,
    )

    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        init_on_cpu = true,
        Courant_number = CFLmax,
    )
    dgn_config = config_diagnostics(driver_config, timeend)

    cbtmarfilter = GenericCallbacks.EveryXSimulationSteps(1) do
        Filters.apply!(
            solver_config.Q,
            ("moisture.ρq_tot",),
            solver_config.dg.grid,
            TMARFilter(),
        )
        nothing
    end

    check_cons = (
        ClimateMachine.ConservationCheck("ρ", "3000steps", FT(0.0001)),
        ClimateMachine.ConservationCheck("energy.ρe", "3000steps", FT(0.0025)),
    )

    result = ClimateMachine.invoke!(
        solver_config;
        user_callbacks = (cbtmarfilter,),
        diagnostics_config = dgn_config,
        check_cons = check_cons,
        check_euclidean_distance = true,
    )
end

function config_diagnostics(driver_config, timeend)
    FT = eltype(driver_config.grid)
    info = driver_config.config_info
    interval = "$(cld(timeend, 2) + 10)ssecs"

    boundaries = [
        FT(0) FT(0) FT(0)
        FT(info.hmax) FT(info.hmax) FT(info.zmax)
    ]
    axes = (
        [FT(1)],
        [FT(1)],
        collect(range(boundaries[1, 3], boundaries[2, 3], step = FT(50)),),
    )
    interpol = ClimateMachine.InterpolationConfiguration(
        driver_config,
        boundaries;
        axes = axes,
    )
    dgngrp = setup_dump_state_diagnostics(
        SingleStackConfigType(),
        interval,
        driver_config.name,
        interpol = interpol,
    )
    return ClimateMachine.DiagnosticsConfiguration([dgngrp])
end

main()
