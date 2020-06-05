#!/usr/bin/env julia --project
using ClimateMachine
ClimateMachine.init()
using ClimateMachine.Atmos
using ClimateMachine.ConfigTypes
using ClimateMachine.Diagnostics
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.SystemSolvers: ManyColumnLU
using ClimateMachine.Mesh.Filters
using ClimateMachine.Mesh.Grids
using ClimateMachine.TemperatureProfiles
using ClimateMachine.Thermodynamics:
    air_temperature, internal_energy, air_pressure
using ClimateMachine.VariableTemplates
using ClimateMachine.MPIStateArrays: realview
using ClimateMachine.ODESolvers: dostep!

using Distributions: Uniform
using LinearAlgebra
using StaticArrays
using Random: rand
using Test

using CLIMAParameters
using CLIMAParameters.Planet: R_d, day, grav, cp_d, cv_d, planet_radius
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

struct HeldSuarezDataConfig{FT}
    T_ref::FT
end

function init_heldsuarez!(bl, state, aux, coords, t)
    FT = eltype(state)

    # Set initial state to reference state with random perturbation
    rnd = FT(1.0 + rand(Uniform(-1e-3, 1e-3)))
    state.ρ = aux.ref_state.ρ
    state.ρu = SVector{3, FT}(0, 0, 0)
    state.ρe = rnd * aux.ref_state.ρe

    nothing
end

function config_heldsuarez(FT, poly_order, resolution, split_explicit_implicit)
    # Set up a reference state for linearization of equations
    temp_profile_ref = DecayingTemperatureProfile{FT}(param_set)
    ref_state = HydrostaticState(temp_profile_ref)

    # Set up a Rayleigh sponge to dampen flow at the top of the domain
    domain_height::FT = 30e3               # distance between surface and top of atmosphere (m)
    z_sponge::FT = 12e3                    # height at which sponge begins (m)
    α_relax::FT = 1 / 60 / 15              # sponge relaxation rate (1/s)
    exp_sponge = 2                         # sponge exponent for squared-sinusoid profile
    u_relax = SVector(FT(0), FT(0), FT(0)) # relaxation velocity (m/s)
    sponge = RayleighSponge{FT}(
        domain_height,
        z_sponge,
        α_relax,
        u_relax,
        exp_sponge,
    )

    # Set up the atmosphere model
    exp_name = "HeldSuarez"
    T_ref::FT = 255        # reference temperature for Held-Suarez forcing (K)
    τ_hyper::FT = 4 * 3600 # hyperdiffusion time scale in (s)
    c_smag::FT = 0.21      # Smagorinsky coefficient
    model = AtmosModel{FT}(
        AtmosGCMConfigType,
        param_set;
        ref_state = ref_state,
        turbulence = SmagorinskyLilly(c_smag),
        hyperdiffusion = StandardHyperDiffusion(τ_hyper),
        moisture = DryModel(),
        source = (Gravity(), Coriolis(), held_suarez_forcing!, sponge),
        init_state_conservative = init_heldsuarez!,
        data_config = HeldSuarezDataConfig(T_ref),
    )

    config = ClimateMachine.AtmosGCMConfiguration(
        exp_name,
        poly_order,
        resolution,
        domain_height,
        param_set,
        init_heldsuarez!;
        model = model,
        solver_type = ClimateMachine.IMEXSolverType(;
            split_explicit_implicit = split_explicit_implicit,
            discrete_splitting = true,
            variant = NaiveVariant(),
        ),
    )

    return config
end

function held_suarez_forcing!(
    bl,
    source,
    state,
    diffusive,
    aux,
    t::Real,
    direction,
)
    FT = eltype(state)

    # Parameters
    T_ref = bl.data_config.T_ref

    # Extract the state
    ρ = state.ρ
    ρu = state.ρu
    ρe = state.ρe

    coord = aux.coord
    e_int = internal_energy(bl.moisture, bl.orientation, state, aux)
    T = air_temperature(bl.param_set, e_int)
    _R_d = FT(R_d(bl.param_set))
    _day = FT(day(bl.param_set))
    _grav = FT(grav(bl.param_set))
    _cp_d = FT(cp_d(bl.param_set))
    _cv_d = FT(cv_d(bl.param_set))

    # Held-Suarez parameters
    k_a = FT(1 / (40 * _day))
    k_f = FT(1 / _day)
    k_s = FT(1 / (4 * _day))
    ΔT_y = FT(60)
    Δθ_z = FT(10)
    T_equator = FT(315)
    T_min = FT(200)
    σ_b = FT(7 / 10)

    # Held-Suarez forcing
    φ = latitude(bl.orientation, aux)
    p = air_pressure(bl.param_set, T, ρ)

    #TODO: replace _p0 with dynamic surfce pressure in Δσ calculations to account
    #for topography, but leave unchanged for calculations of σ involved in T_equil
    _p0 = 1.01325e5
    σ = p / _p0
    exner_p = σ^(_R_d / _cp_d)
    Δσ = (σ - σ_b) / (1 - σ_b)
    height_factor = max(0, Δσ)
    T_equil = (T_equator - ΔT_y * sin(φ)^2 - Δθ_z * log(σ) * cos(φ)^2) * exner_p
    T_equil = max(T_min, T_equil)
    k_T = k_a + (k_s - k_a) * height_factor * cos(φ)^4
    k_v = k_f * height_factor

    # Apply Held-Suarez forcing
    source.ρu -= k_v * projection_tangential(bl, aux, ρu)
    source.ρe -= k_T * ρ * _cv_d * (T - T_equil)
    return nothing
end

function config_diagnostics(FT, driver_config)
    interval = "40000steps" # chosen to allow a single diagnostics collection

    _planet_radius = FT(planet_radius(param_set))

    info = driver_config.config_info
    boundaries = [
        FT(-90.0) FT(-180.0) _planet_radius
        FT(90.0) FT(180.0) FT(_planet_radius + info.domain_height)
    ]
    resolution = (FT(10), FT(10), FT(1000)) # in (deg, deg, m)
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

function main()
    # Driver configuration parameters
    FT = Float64                             # floating type precision
    poly_order = 5                           # discontinuous Galerkin polynomial order
    n_horz = 5                               # horizontal element number
    n_vert = 5                               # vertical element number
    n_days = 50                              # experiment day number
    timestart = FT(0)                        # start time (s)
    timeend = FT(n_days * day(param_set))    # end time (s)

    # Set up driver configuration
    driver_config_false =
        config_heldsuarez(FT, poly_order, (n_horz, n_vert), false)
    driver_config_true =
        config_heldsuarez(FT, poly_order, (n_horz, n_vert), true)

    # Set up experiment
    solver_config_false = ClimateMachine.SolverConfiguration(
        timestart,
        timeend,
        driver_config_false,
        Courant_number = 0.2,
        init_on_cpu = true,
        CFL_direction = HorizontalDirection(),
        diffdir = HorizontalDirection(),
    )

    # Set up experiment
    solver_config_true = ClimateMachine.SolverConfiguration(
        timestart,
        timeend,
        driver_config_true,
        Courant_number = 0.2,
        init_on_cpu = true,
        CFL_direction = HorizontalDirection(),
        diffdir = HorizontalDirection(),
    )

    Q = solver_config_false.Q

    Ql = Array(Q.data)
    dostep!(Q, solver_config_false.solver, nothing, FT(0))
    Qfalse = Array(realview(Q))
    copy!(Q.data, Ql)
    dostep!(Q, solver_config_true.solver, nothing, FT(0))
    Qtrue = Array(realview(Q))
    dQ = (Qfalse - Qtrue)
    normalized_dQ = dQ ./ max.(1, max.(abs.(Qfalse), abs.(Qtrue)))
    println()
    @show extrema(dQ)
    @show extrema(normalized_dQ)

    QS_true = solver_config_true.solver.Qstages
    QS_false = solver_config_false.solver.Qstages
    for s in 1:length(QS_true)
        QT = Array(QS_true[s])
        QF = Array(QS_false[s])
        dQ = (QT - QF)
        normalized_dQ = dQ ./ max.(1, max.(abs.(QT), abs.(QF)))
        println()
        @show s
        @show extrema(dQ)
        @show extrema(normalized_dQ)
    end

    Ql .+= (2 * rand(size(Ql)...) .- 1) .* Ql / 5
    copy!(Q.data, Ql)

    full = solver_config_false.solver.rhs!
    full2 = solver_config_true.solver.rhs!

    dQfull = similar(Q)
    dQfull2 = similar(Q)

    full(dQfull, Q, nothing, FT(0), increment = false)
    full2(dQfull2, Q, nothing, FT(0), increment = false)

    Afull = Array(realview(dQfull))
    Afull2 = Array(realview(dQfull2))

    dA = (Afull - Afull2)
    normalized_dA = dA ./ max.(1, max.(abs.(Afull), abs.(Afull2)))
    println()
    @show extrema(dA)
    @show extrema(normalized_dA)

end

main()
nothing
