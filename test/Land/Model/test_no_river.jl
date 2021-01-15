# Test that the land model still runs, even with the lowest/simplest
# version of soil (prescribed heat and prescribed water - no state
# variables)
using MPI
using OrderedCollections
using StaticArrays
using Test
using Statistics

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using ClimateMachine
using ClimateMachine.Land
using ClimateMachine.Land.River
using ClimateMachine.Land.SoilWaterParameterizations
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

@testset "NoRiver Model" begin
    ClimateMachine.init()
    FT = Float64

    function init_land_model!(land, state, aux, localgeo, time) end

    soil_water_model = PrescribedWaterModel()
    soil_heat_model = PrescribedTemperatureModel()
    soil_param_functions = nothing

    m_soil = SoilModel(soil_param_functions, soil_water_model, soil_heat_model)
    m_river = NoRiverModel()

    sources = ()
    m = LandModel(
        param_set,
        m_soil,
        m_river;
        source = sources,
        init_state_prognostic = init_land_model!,
    )

    N_poly = 5
    nelem_vert = 10

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
    timeend = FT(60)
    dt = FT(1)

    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        ode_dt = dt,
    )
    mygrid = solver_config.dg.grid
    Q = solver_config.Q
    aux = solver_config.dg.state_auxiliary

    ClimateMachine.invoke!(solver_config;)

    t = ODESolvers.gettime(solver_config.solver)
    state_vars = SingleStackUtils.get_vars_from_nodal_stack(
        mygrid,
        Q,
        vars_state(m, Prognostic(), FT),
    )
    aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
        mygrid,
        aux,
        vars_state(m, Auxiliary(), FT),
    )
    #Make sure it runs, and that there are no state variables, and only "z" as aux.
    @test t == timeend
    @test size(Base.collect(keys(aux_vars)))[1] == 3
    @test size(Base.collect(keys(state_vars)))[1] == 0
end


@testset "Analytical River Model" begin
    ClimateMachine.init()
    FT = Float64

    function init_land_model!(land, state, aux, localgeo, time) end

    soil_water_model = PrescribedWaterModel()
    soil_heat_model = PrescribedTemperatureModel()
    soil_param_functions = nothing

    m_soil = SoilModel(soil_param_functions, soil_water_model, soil_heat_model)
    m_river = RiverModel(
        (x,y) -> eltype(x)(1),
        (x,y) -> eltype(x)(0.0),
        (x,y) -> eltype(x)(0.0016),
        (x,y) -> eltype(x)(1);
        mannings = (x,y) -> eltype(x)(0.025)
    )

    # 10.1061/(ASCE)0733-9429(2007)133:2(217) 
    # Eqn 6
    # debug boundry condition, constant positive flow rate after min 30
    # add river boundry state land domain bc
    bc = LandDomainBC(
        minx_bc = LandComponentBC(river = Dirichlet((aux, t) -> eltype(aux)(0))),
    )
 
    function init_land_model!(land, state, aux, localgeo, time)
        state.river.area = eltype(state)(0)
    end

    # units in m / s 
    precip(x, y, t) = t < (30 * 60) ? 1.4e-5 : 0.0

    sources = (Precip{FT}(precip),)

    m = LandModel(
        param_set,
        m_soil,
        m_river;
        boundary_conditions = bc,
        source = sources,
        init_state_prognostic = init_land_model!,
    )
    function warp_constant_slope(xin, yin, zin; topo_max = 0.2, zmin = -0.1, xmax = 400)
        FT = eltype(xin)
        zmax = FT((FT(1.0)-xin / xmax) * topo_max)
        alpha = FT(1.0) - zmax / zmin
        zout = zmin + (zin - zmin) * alpha
        x, y, z = xin, yin, zout
        return x, y, z
    end

    N_poly = 1;
    xres = FT(2.286)
    yres = FT(0.25)
    zres = FT(0.1)
    # Specify the domain boundaries.
    zmax = FT(0);
    zmin = FT(-0.1);
    xmax = FT(182.88)
    ymax = FT(1.0)
    topo_max = FT(0.0016*xmax)

    driver_config = ClimateMachine.MultiColumnLandModel(
        "LandModel",
        (N_poly, N_poly),
        (xres,yres,zres),
        xmax,
        ymax,
        zmax,
        param_set,
        m;
        zmin = zmin,
        #numerical_flux_first_order = CentralNumericalFluxFirstOrder(),now the default for us
        meshwarp = (x...) -> warp_constant_slope(x...;
        topo_max = topo_max, zmin = zmin, xmax = xmax),
    );

    t0 = FT(0)
    timeend = FT(60*60)
    dt = FT(1)

    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        ode_dt = dt,
    )
    mygrid = solver_config.dg.grid
    Q = solver_config.Q
    
    area_index =
        varsindex(vars_state(m, Prognostic(), FT), :river, :area)
    n_outputs = 60

    every_x_simulation_time = ceil(Int, timeend / n_outputs)

    dons = Dict([k => Dict() for k in 1:n_outputs]...)

    iostep = [1]
    callback = GenericCallbacks.EveryXSimulationTime(
        every_x_simulation_time,
    ) do (init = false)
        t = ODESolvers.gettime(solver_config.solver)
        area = Q[:, area_index, :]
        all_vars = Dict{String, Array}(
            "t" => [t],
            "area" => area,
        )
        dons[iostep[1]] = all_vars
        iostep[1] += 1
        return
    end

    ClimateMachine.invoke!(solver_config; user_callbacks = (callback,))


    aux = solver_config.dg.state_auxiliary;
    mask = aux[:,1,:] .== 182.88
    n_outputs = length(dons)
    area = [mean(dons[k]["area"][mask[:]]) for k in 1:n_outputs]
    height = area ./ ymax
    time_data = [dons[l]["t"][1] for l in 1:n_outputs]
    
    alpha = sqrt(0.0016)/0.025
    i = 1.4e-5
    L = xmax
    m = 5/3
    t_c = (L*i^(1-m)/alpha)^(1/m)
    t_r = 30*60
    
    q = height.^(m) .* alpha
    
    function g(m,y, i, t_r, L, alpha, t)
        output = L/alpha-y^(m)/i-y^(m-1)*m*(t-t_r)
        return output
    end
    
    function dg(m,y, i, t_r, L, alpha, t)
        output = -y^(m-1)*m/i-y^(m-2)*m*(m-1)*(t-t_r)
        return output
    end
    
    function analytic(t,alpha, t_c, t_r, i, L, m)
        if t < t_c
            return alpha*(i*t)^(m)
        end
        
        if t <= t_r && t > t_c
            return alpha*(i*t_c)^(m)
        end
        
        if t > t_r
            yL = (i*(t-t_r))
            delta = 1
            error = g(m,yL,i,t_r,L,alpha,t)
            while abs(error) > 1e-4
                delta = -g(m,yL,i,t_r,L,alpha,t)/dg(m,yL,i,t_r,L,alpha,t)
                yL = yL+ delta
                error = g(m,yL,i,t_r,L,alpha,t)
            end
            return alpha*yL^m    
            
        end
        
    end
    
    @test sqrt_rmse_over_max_q = sqrt(mean((analytic.(time_data, Ref(alpha), Ref(t_c), Ref(t_r), Ref(i), Ref(L), Ref(m)) .- q).^2.0))/ maximum(q) < 3e-3
    
end
