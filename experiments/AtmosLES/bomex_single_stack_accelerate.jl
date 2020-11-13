include("bomex_model.jl")
using ClimateMachine.SingleStackUtils
using ClimateMachine.SystemSolvers
using ClimateMachine.Atmos: vars_state
const clima_dir = dirname(dirname(pathof(ClimateMachine)));
include(joinpath(clima_dir, "docs", "plothelpers.jl"))
ENV["CLIMATEMACHINE_SETTINGS_MONITOR_COURANT_NUMBERS"] = "3000steps"
ENV["CLIMATEMACHINE_SETTINGS_MONITOR_TIMESTEP_DURATION"] = "3000steps"
ENV["CLIMATEMACHINE_SETTINGS_FIX_RNG_SEED"] = true

function main(::Type{FT}) where {FT}
    # add a command line argument to specify the kind of surface flux
    # TODO: this will move to the future namelist functionality
    bomex_args = ArgParseSettings(autofix_names = true)
    add_arg_group!(bomex_args, "BOMEX")
    @add_arg_table! bomex_args begin
        "--surface-flux"
        help = "specify surface flux for energy and moisture"
        metavar = "prescribed|bulk"
        arg_type = String
        default = "prescribed"
    end

    cl_args =
        ClimateMachine.init(parse_clargs = true, custom_clargs = bomex_args)

    surface_flux = cl_args["surface_flux"]

    config_type = SingleStackConfigType

    # DG polynomial order
    N = 3

    # Prescribe domain parameters
    nelem_vert = 15
    zmax = FT(3000)

    t0 = FT(0)

    # For a full-run, please set the timeend to 3600*6 seconds
    # For the test we set this to == 30 minutes
    timeend = FT(1800 * 2)
    #timeend = FT(3600 * 6)
    CFLmax = FT(2)

    # Choose default IMEX solver
    # ode_solver_type = ClimateMachine.IMEXSolverType()

    ode_solver_type = ClimateMachine.IMEXSolverType(
        implicit_model = AtmosAcousticGravityLinearModel,
        implicit_solver = SingleColumnLU,
        solver_method = ARK2GiraldoKellyConstantinescu,
        split_explicit_implicit = true,
        discrete_splitting = false,
    )


    model = bomex_model(FT, config_type, zmax, surface_flux)
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
    dgn_config = config_diagnostics(driver_config)

    state_types = (Prognostic(), Auxiliary())
    dons_arr = [dict_of_nodal_states(solver_config, state_types; interp = true)]
    time_data = FT[0]

    # Define the number of outputs from `t0` to `timeend`
    n_outputs = 10
    # This equates to exports every ceil(Int, timeend/n_outputs) time-step:
    every_x_simulation_time = ceil(Int, timeend / n_outputs)

    cb_data_vs_time =
        GenericCallbacks.EveryXSimulationTime(every_x_simulation_time) do
            push!(dons_arr, dict_of_nodal_states(solver_config, state_types; interp = true))
            push!(time_data, gettime(solver_config.solver))
            nothing
        end

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
        ClimateMachine.ConservationCheck("ρe", "3000steps", FT(0.0025)),
    )

    result = ClimateMachine.invoke!(
        solver_config;
        user_callbacks = (cbtmarfilter, cb_data_vs_time),
        diagnostics_config = dgn_config,
        check_cons = check_cons,
        check_euclidean_distance = true,
    )

    dons = dict_of_nodal_states(solver_config, state_types; interp = true)
    push!(dons_arr, dons)
    push!(time_data, gettime(solver_config.solver))
    return solver_config, dons_arr, time_data, state_types
end

solver_config, dons_arr, time_data, state_types = main(Float64)

export_state_plots(
    solver_config,
    dons_arr,
    time_data,
    joinpath("output", "bomex_ss_acc");
    z = Array(get_z(solver_config.dg.grid; rm_dupes = true)),
)
