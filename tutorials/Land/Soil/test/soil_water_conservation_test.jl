# --------------------------------- run CLIMA SOIL MODEL -----------------------
# CLIMA_run_SoilWater.jl: This model simulates soil water dynamics for the CliMA model

println("1) Import/Export Needed Functions")

# Load necessary CliMA subroutines
using MPI
using OrderedCollections
using Plots
using StaticArrays
using CLIMAParameters
using DocStringExtensions
using CLIMAParameters.Planet: day
using CLIMAParameters.Planet: ρ_cloud_liq
struct EarthParameterSet <: AbstractEarthParameterSet end
const Earth_param_set = EarthParameterSet()
# struct LandParameterSet <: AbstractLandParameterSet end
# const Land_param_set = LandParameterSet()
using ClimateMachine
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.Writers
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.DGMethods: BalanceLaw, LocalGeometry
using ClimateMachine.MPIStateArrays
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.VariableTemplates
using ClimateMachine.SingleStackUtils
using DelimitedFiles
using CSV
using Parameters
#we need to add these since we want to do integrals
import ClimateMachine.DGMethods:
    vars_integrals,
    integral_load_auxiliary_state!,
    integral_set_auxiliary_state!,
    indefinite_stack_integral!
include("soil_model_base_numerics.jl")
include("conservation_numerics.jl")
# Add water functions
include("../Water/water_functions.jl")

FT = Float64

# Initialize CliMA
ClimateMachine.init(; disable_gpu = true);

# Determine formulation to use for K and ψ_m
WF = waterfunctions(
    hydraulic_cond = Haverkamp{FT}(),
    matric_pot = vanGenuchten{FT}()
)

# Constants from Parameters Module
_day = FT(day(Earth_param_set))
_hour = _day/24
_minute  = _hour/60
_ρ_cloud_liq = FT(ρ_cloud_liq(Earth_param_set))

# Constants which will change in space, to import from database
porosity = 0.495 # m3/m3 Read in from data base
S_s = 10e-4  # m-1 Read in from data base
K_sat  = 0.0443 / (3600*100) # "Clay" Read in from data base?
#K_sat  = 34 / (3600*100) # "Sand"

println("2) Prep ICs, and time-stepper and output configurations...")
# Configure a `ClimateMachine` solver.
const timeend = FT(_hour)
const t0 = FT(0)

println("3) Set up system...")
#IC/BC values
ν_0 = 0.24
ν_surface = porosity - 1e-3
ψ_0 = pressure_head(WF.matric_pot,porosity,S_s,ν_0)
S_l_0 = effective_saturation(porosity, ν_0)
κ_0 = hydraulic_conductivity(WF.hydraulic_cond, K_sat, S_l_0, ψ_0)
#Note - this is actualy grad h, not flux (grad h *kappa)
surface_flux = 0.5
bottom_flux = -1.0
# Load Soil Model in 'm'
m = SoilModelMoisture(
    param_set = Earth_param_set,
    WF = WF,
    initialκ   = (aux) -> κ_0,
    initialν = (state, aux) -> ν_0,
    initialh = (aux) -> aux.z + ψ_0,
    dirichlet_bc = Dirichlet(NaN, NaN),
    neumann_bc = Neumann(surface_flux, bottom_flux)
)

println("4) Define variables for simulation...") # move up

# # Spatial discretization

# Prescribe polynomial order of basis functions in finite elements
N_poly = 5;

# Specify the number of vertical elements
nelem_vert = 10;

# Specify the domain boundaries
zmax = FT(0);
zmin = FT(-1)

# Establish a `ClimateMachine` single stack configuration
driver_config = ClimateMachine.SingleStackConfiguration(
    "SoilMoistureModel",
    N_poly,
    nelem_vert,
    zmax,
    Earth_param_set,
    m;
    zmin = zmin,
    numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
);

# Minimum spatial and temporal steps
Δ = min_node_distance(driver_config.grid)
τ = (Δ^2 /K_sat)# this should be Kmax
dt = 6# TBD - make automatic^

# This initializes the state vector and allocates memory for the solution in
# space (`dg` has the model `m`, which describes the PDEs as well as the
# function used for initialization). This additionally initializes the ODE
# solver, by default an explicit Low-Storage
# [Runge-Kutta](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
# method.

solver_config =
    ClimateMachine.SolverConfiguration(t0, timeend, driver_config, ode_dt = dt);
mygrid = solver_config.dg.grid;
Q = solver_config.Q;
aux = solver_config.dg.state_auxiliary;

# # Solver hooks / callbacks

# Define the number of outputs from `t0` to `timeend`. Note that t = 0 is considered an output, so we will need to get the output after the integration is done as well.
const n_outputs = 10;

# This equates to exports every ceil(Int, timeend/n_outputs) time-step:
const every_x_simulation_time = ceil(Int, timeend / n_outputs);

# Create a nested dictionary to store the solution:
all_data = Dict([k => Dict() for k in 0:n_outputs]...)
# The `ClimateMachine`'s time-steppers provide hooks, or callbacks, which
# allow users to inject code to be executed at specified intervals.
step = [0];
callback = GenericCallbacks.EveryXSimulationTime(
    every_x_simulation_time,
    solver_config.solver,
) do (init = false)
    t = ODESolvers.gettime(
        solver_config.solver
    )
    grads = SingleStackUtils.get_vars_from_nodal_stack(
        mygrid,
        solver_config.dg.state_gradient_flux,
        vars_state_gradient_flux(m, FT),
    )    
    state_vars = SingleStackUtils.get_vars_from_nodal_stack(
        mygrid,
        Q,
        vars_state_conservative(m, FT),
    )
    aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
        mygrid,
        aux,
        vars_state_auxiliary(m, FT);
        #exclude = ["z"],
    )
    all_vars = OrderedDict(state_vars..., aux_vars..., grads...)
    all_vars["t"]= [t]
    all_data[step[1]] = all_vars

    step[1] += 1
    nothing
end;

# # Solve

ClimateMachine.invoke!(solver_config; user_callbacks = (callback,));
#Pull out state etc. at final time step
t = ODESolvers.gettime(
    solver_config.solver
)
grads = SingleStackUtils.get_vars_from_nodal_stack(
    mygrid,
    solver_config.dg.state_gradient_flux,
    vars_state_gradient_flux(m, FT),
)
state_vars = SingleStackUtils.get_vars_from_nodal_stack(
    mygrid,
    Q,
    vars_state_conservative(m, FT),
)
aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
    mygrid,
    aux,
    vars_state_auxiliary(m, FT);
)
all_vars = OrderedDict(state_vars..., aux_vars..., grads...);
all_vars["t"]= [t]
all_data[n_outputs] = all_vars


#Compute expected and actual energy changes as a function of time, using prescribed flux and source terms.
mass_total = [all_data[i]["int.a"][end] for i in 0:1:n_outputs]
mass_total = mass_total.*_ρ_cloud_liq

times = [all_data[i]["t"][end] for i in 0:1:n_outputs]
hydraulic_k_bottom = [all_data[i]["κ"][1] for i in 0:1:n_outputs]
hydraulic_k_top = [all_data[i]["κ"][end] for i in 0:1:n_outputs]

#check signs on normal vectors in flux pieces.
integrated_surface_flux = surface_flux*1.0.*hydraulic_k_top.*times
integrated_bottom_flux = bottom_flux*(-1.0).*hydraulic_k_bottom.*times

net_integrated_flux = (integrated_surface_flux.-integrated_bottom_flux).*_ρ_cloud_liq
domain_volume = 1.0
initial_mass = ν_0*domain_volume*_ρ_cloud_liq
analytic_mass = (initial_mass.+net_integrated_flux)
y = (mass_total.-analytic_mass)
worst_error = maximum(abs.(y))

#assert worst_error is < some threshold TBD
#add in sources
#allow temporal dependency in BC?

#plot(times, log10.(abs.(mass_total.-initial_mass)),seriestype = :scatter, markserize = 8, xlabel = "Time (s)", ylabel = "log10(|Δm|)", label = "numerical", ylims = [-8,-4])
#plot!(times, log10.(abs.(analytic_mass.-initial_mass)),seriestype = :scatter, markersize = 2,label = "expected")

#plot(times, log10.(abs.(mass_total.-analytic_mass)),seriestype = :scatter, markserize = 8, xlabel = "Time (s)", ylabel = "log10(|Δm|)", label = "numerical")
