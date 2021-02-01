# Test that the way we specify boundary conditions works as expected
using MPI
using OrderedCollections
using StaticArrays
using Statistics
using Test

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using ClimateMachine
using ClimateMachine.Land
using ClimateMachine.Land.SoilWaterParameterizations
using ClimateMachine.Land.Runoff
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.DGMethods: BalanceLaw, LocalGeometry
using ClimateMachine.MPIStateArrays
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.VariableTemplates
using ClimateMachine.SingleStackUtils
using ClimateMachine.BalanceLaws:
    BalanceLaw, Prognostic, Auxiliary, Gradient, GradientFlux, vars_state

import ClimateMachine.DGMethods.FVReconstructions: FVLinear

@testset "NoRunoff" begin
    ClimateMachine.init()

    FT = Float64

    function init_soil_water!(land, state, aux, localgeo, time)
        myfloat = eltype(state)
        state.soil.water.ϑ_l = myfloat(land.soil.water.initialϑ_l(aux))
        state.soil.water.θ_i = myfloat(land.soil.water.initialθ_i(aux))
    end

    soil_param_functions =
        SoilParamFunctions{FT}(porosity = 0.75, Ksat = 1e-7, S_s = 1e-3)
    surface_precip_amplitude = FT(3e-8)
    f = FT(pi * 2.0 / 300.0)
    precip = (t) -> surface_precip_amplitude * sin(f * t)
    ϑ_l0 = (aux) -> eltype(aux)(0.2)
    bc = LandDomainBC(
        bottom_bc = LandComponentBC(soil_water = Neumann((aux,t) -> eltype(aux)(0.0))),
        surface_bc = LandComponentBC( soil_water = SurfaceDrivenWaterBoundaryConditions(
            FT;
            precip_model = DrivenConstantPrecip{FT}(precip),
            runoff_model = NoRunoff(),
        ))
    )

    soil_water_model = SoilWaterModel(FT; initialϑ_l = ϑ_l0)
    soil_heat_model = PrescribedTemperatureModel()

    m_soil = SoilModel(soil_param_functions, soil_water_model, soil_heat_model)
    sources = ()
    m = LandModel(
        param_set,
        m_soil;
        boundary_conditions = bc,
        source = sources,
        init_state_prognostic = init_soil_water!,
    )


    N_poly = 5
    nelem_vert = 50


    # Specify the domain boundaries
    zmax = FT(0)
    zmin = FT(-1)

    driver_config = ClimateMachine.SingleStackConfiguration(
        "LandModel",
        N_poly,
        nelem_vert,
        zmax,
        param_set,
        m;
        zmin = zmin,
        numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
    )


    t0 = FT(0)
    timeend = FT(150)
    dt = FT(0.05)

    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        ode_dt = dt,
    )
    n_outputs = 60
    mygrid = solver_config.dg.grid
    const every_x_simulation_time = ceil(Int, timeend / n_outputs);
    state_types = (Prognostic(), Auxiliary(), GradientFlux())
    dons_arr = Dict[dict_of_nodal_states(solver_config, state_types; interp = true)]
    time_data = FT[0] # store time data
    
    # We specify a function which evaluates `every_x_simulation_time` and returns
    # the state vector, appending the variables we are interested in into
    # `dons_arr`.
    
    callback = GenericCallbacks.EveryXSimulationTime(every_x_simulation_time) do
        dons = dict_of_nodal_states(solver_config, state_types; interp = true)
        push!(dons_arr, dons)
        push!(time_data, gettime(solver_config.solver))
        nothing
    end;
    
    # # Run the integration
    ClimateMachine.invoke!(solver_config; user_callbacks = (callback,));

    dons = dict_of_nodal_states(solver_config, state_types; interp = true)
    push!(dons_arr, dons)
    push!(time_data, gettime(solver_config.solver));
    
    # Get z-coordinate
    z = get_z(solver_config.dg.grid; rm_dupes = true);
    N = length(dons_arr)
    computed_surface_flux = [dons_arr[k]["soil.water.K∇h[3]"][end] for k in 1:N]
    t = time_data
    prescribed_surface_flux = t -> FT(-1) * FT(3e-8 * sin(pi * 2.0 * t / 300.0))
    MSE = mean((prescribed_surface_flux.(t) .- computed_surface_flux) .^ 2.0)
    @test MSE < 1e-7
end


#=
@testset "Infiltration Excess" begin
    ClimateMachine.init()

    FT = Float64

    function init_soil_water!(land, state, aux, localgeo, time)
        myfloat = eltype(state)
        state.soil.water.ϑ_l = myfloat(land.soil.water.initialϑ_l(aux))
        state.soil.water.θ_i = myfloat(land.soil.water.initialθ_i(aux))
    end

    soil_param_functions =
        SoilParamFunctions{FT}(porosity = 0.495, Ksat = 0.0443 / (3600 * 100), S_s = 1e-4)
    surface_precip_amplitude = FT(-5e-7)
    f = FT(pi * 2.0 / 2000.0)
    function precip_function(t, f)
        myf = eltype(t)
        if t < myf(10000.0)
            return myf(0.0)
        elseif t < myf(10500.0)
            return surface_precip_amplitude * sin(f * (t-myf(10000)))
        else
            return surface_precip_amplitude
        end
    end

    
    precip = (t) -> precip_function(t, f)
    hydraulics = vanGenuchten{FT}(n = 2.0)
    function hydrostatic_profile(z, zm, porosity, n, α)
        myf = eltype(z)
        m = FT(1-1/n)
        S = FT((FT(1)+(α*(z-zm))^n)^(-m))
        return FT(S*porosity)
    end

    


    ϑ_l0 = (aux) -> eltype(aux)(hydrostatic_profile(aux.z,-1,0.495, hydraulics.n, hydraulics.α))
    bc = LandDomainBC(
        bottom_bc = LandComponentBC(soil_water = Neumann((aux,t) -> eltype(aux)(0.0))),
        surface_bc = LandComponentBC( soil_water = SurfaceDrivenWaterBoundaryConditions(
            FT;
            precip_model = DrivenConstantPrecip{FT}(precip),
            runoff_model = CoarseGridRunoff(),
        ))
    )

    soil_water_model = SoilWaterModel(
        FT;
        moisture_factor = MoistureDependent{FT}(),
        hydraulics = hydraulics,
        initialϑ_l = ϑ_l0,
    );
    soil_heat_model = PrescribedTemperatureModel()

    m_soil = SoilModel(soil_param_functions, soil_water_model, soil_heat_model)
    sources = ()
    m = LandModel(
        param_set,
        m_soil;
        boundary_conditions = bc,
        source = sources,
        init_state_prognostic = init_soil_water!,
    )

    N_poly = (1,0)
    nelem_vert = 100

    # Specify the domain boundaries
    zmax = FT(0)
    zmin = FT(-1)

    driver_config = ClimateMachine.SingleStackConfiguration(
        "LandModel",
        N_poly,
        nelem_vert,
        zmax,
        param_set,
        m;
        zmin = zmin,
        numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
        fv_reconstruction = FVLinear(),

    )


    t0 = FT(0)
    timeend = FT(64000)
    dt = FT(2)

    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        ode_dt = dt,
    )
    n_outputs = 60
    mygrid = solver_config.dg.grid
    const every_x_simulation_time = ceil(Int, timeend / n_outputs);
    state_types = (Prognostic(), Auxiliary(), GradientFlux())
    dons_arr = Dict[dict_of_nodal_states(solver_config, state_types; interp = true)]
    time_data = FT[0] # store time data
    
    # We specify a function which evaluates `every_x_simulation_time` and returns
    # the state vector, appending the variables we are interested in into
    # `dons_arr`.
    
    callback = GenericCallbacks.EveryXSimulationTime(every_x_simulation_time) do
        dons = dict_of_nodal_states(solver_config, state_types; interp = true)
        push!(dons_arr, dons)
        push!(time_data, gettime(solver_config.solver))
        nothing
    end;
    
    # # Run the integration
    ClimateMachine.invoke!(solver_config; user_callbacks = (callback,));

    dons = dict_of_nodal_states(solver_config, state_types; interp = true)
    push!(dons_arr, dons)
    push!(time_data, gettime(solver_config.solver));
    
    # Get z-coordinate
    z = get_z(solver_config.dg.grid; rm_dupes = true);
    N = length(dons_arr)
    computed_surface_flux = [dons_arr[k]["soil.water.K∇h[3]"][end] for k in 1:N]
    t = time_data
end
=#
