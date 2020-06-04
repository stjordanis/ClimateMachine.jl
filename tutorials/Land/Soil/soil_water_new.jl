# --------------------------------- run CLIMA SOIL MODEL -----------------------
# CLIMA_run_SoilWater.jl: This model simulates soil water dynamics for the CliMA model

######
###### 1) Import/Export Needed Functions
######
println("1) Import/Export Needed Functions")

# Load necessary CliMA subroutines
using MPI
using OrderedCollections
using Plots
using StaticArrays
using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const my_param_set = EarthParameterSet()

#using Test
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

#using Logging
#using Printf
#using NCDatasets

#import ClimateMachine.DGMethods:
#    vars_state_auxiliary,
#    vars_state_conservative,
#    vars_state_gradient,
#    vars_state_gradient_flux,
#    source!,
#    flux_second_order!,
#    flux_first_order!,
#    compute_gradient_argument!,
#    compute_gradient_flux!,
#    update_auxiliary_state!,
#    nodal_update_auxiliary_state!,
#    init_state_auxiliary!,
#    init_state_conservative!,
#    boundary_state!

# ENV["GKS_ENCODING"] = "utf-8"
FT = Float64
# Initialize CliMA
ClimateMachine.init(; disable_gpu = true);
const clima_dir = dirname(dirname(pathof(ClimateMachine)));
include(joinpath(clima_dir, "docs", "plothelpers.jl"));


# Change output directory and save plots there
output_dir = joinpath(dirname(dirname(pathof(ClimateMachine))), "output", "land")
mkpath(output_dir)

# Add soil moisture model
include("soil_water_model_new.jl")

# Add water functions
#include("Water/soil_water_properties.jl")
#include("Water/frozen_impedence_factor.jl")
#include("Water/temperature_dependence.jl")
include("Water/matric_potential.jl")
include("Water/pressure_head.jl")
include("Water/hydraulic_head.jl")
include("Water/effective_saturation.jl")
include("Water/hydraulic_conductivity.jl")
#include("Water/augmented_liquid.jl")
#include("Water/calculate_frozen_water.jl")
#include("Water/heaviside.jl")

######
###### 2) Set up system
######
println("2) Set up system...")

# Read in state variables and data
mineral_properties = "Clay"
#soil_T = 280 # Read in from heat model {aux.T}
#soil_Tref = 282.42 # Soil reference temperature: annual mean temperature of site
#theta_liq_0 = 0.2 # Read in from water model {state.θ}
#theta_liq_surface = 0.2 # Read in from water model {state.θ}
#theta_ice_0 = 0 # Read in from water model {state.θi}
#h_0 = -30 # Read in from water model {state.θ}
#ψ_0 = -20 # Soil pressure head {aux.h}
K_sat  = 0.0443 / (3600*100)
porosity = 0.495 # Read in from data base
S_s = 10e-4  # [ m-1] ## Check this value !!
flag = "van Genuchten" # "van Genuchten" , "Brooks and Corey"
ν_0 = 0.24
ν_surface = porosity-1e-3
S_l_0 = effective_saturation(porosity, ν_0)
ψ_0 = pressure_head(S_l_0,porosity,S_s,ν_0,flag)
println(ψ_0)
κ_0 = hydraulic_conductivity(K_sat, S_l_0, ψ_0,0, "Havercamp")# only tested havercamp so far. van Genuchten required very small time step and still had Domain error. Note that this function is called in the solver too, since κ is an auxilary variable. we can think about this, non ideal to define it in both places, also per Elias may not be ideal to have κ be auxiliary.

#This is implemented in hydraulic_conductivity now, with flag "Havercamp".
#function ksat_function(K_sat, head,zed)
#    return K_sat*124.6/(124.6+abs(100*(head-zed))^1.77)
#end


# Load Soil Model in 'm'
m = SoilModelMoisture(
    param_set = my_param_set,
    # Define hydraulic conductivity of soil
    #It seems like initial aux variable values cant depend on state.
    initialκ   = (aux) -> κ_0,
    # Define initial soil moisture
    initialν = (state, aux) -> ν_0, #theta_liq_0, # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287,
    surfaceν = (state, aux, t) -> ν_surface, #theta_liq_surface, # [m3/m3] constant flux at surface, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
    # Define initial and boundary condition parameters
    #initialS_l = (aux) -> S_l_0,
    #initialψ_m = (aux) -> ψ_m_0,
    #initialψ = (aux) -> ψ_0,
    initialh = (aux) -> aux.z + ψ_0 # [m3/m3] constant water content in soil
)

######
###### 3) Define variables for simulation
######
println("3) Define variables for simulation...") # move up


# # Spatial discretization

# Prescribe polynomial order of basis functions in finite elements
N_poly = 5;

# Specify the number of vertical elements
nelem_vert = 10;

# Specify the domain height
zmax = FT(1);


# Establish a `ClimateMachine` single stack configuration
driver_config = ClimateMachine.SingleStackConfiguration(
    "SoilMoistureModel",
    N_poly,
    nelem_vert,
    zmax,
    my_param_set,
    m,
    numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
);


# Minimum spatial and temporal steps
Δ = min_node_distance(driver_config.grid)
τ = (Δ^2 /K_sat)
dt = 0.01*τ #CFL_bound*0.5 # TODO: provide a "default" timestep based on  Δx,Δy,Δz


# Define time variables
const minute = 60
const hour = 60*minute
const day = 24*hour
# const timeend = 1*minute
# const n_outputs = 25
const timeend = FT(2*day)
const t0 = FT(0)

######
###### 4) Prep ICs, and time-stepper and output configurations
######
println("4) Prep ICs, and time-stepper and output configurations...")
# # Configure a `ClimateMachine` solver.

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


state_vars = SingleStackUtils.get_vars_from_nodal_stack(
    mygrid,
    Q,
    vars_state_conservative(m, FT),
)
aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
    mygrid,
    aux,
    vars_state_auxiliary(m, FT),
)
all_vars = OrderedDict(state_vars..., aux_vars...);


# # Solver hooks / callbacks

# Define the number of outputs from `t0` to `timeend`
const n_outputs = 5;

# This equates to exports every ceil(Int, timeend/n_outputs) time-step:
const every_x_simulation_time = ceil(Int, timeend / n_outputs);

# Create a nested dictionary to store the solution:
all_data = Dict([k => Dict() for k in 0:n_outputs]...)
all_data[0] = all_vars # store initial condition at ``t=0``

# The `ClimateMachine`'s time-steppers provide hooks, or callbacks, which
# allow users to inject code to be executed at specified intervals. In this
# callback, the state and aux variables are collected, combined into a single
# `OrderedDict` and written to a NetCDF file (for each output step `step`).
step = [1];
callback = GenericCallbacks.EveryXSimulationTime(
    every_x_simulation_time,
    solver_config.solver,
) do (init = false)
    state_vars = SingleStackUtils.get_vars_from_nodal_stack(
        mygrid,
        Q,
        vars_state_conservative(m, FT),
    )
    aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
        mygrid,
        aux,
        vars_state_auxiliary(m, FT);
        exclude = ["z"],
    )
    all_vars = OrderedDict(state_vars..., aux_vars...)
    all_data[step[1]] = all_vars

    step[1] += 1
    nothing
end;

# # Solve

# This is the main `ClimateMachine` solver invocation. While users do not have
# access to the time-stepping loop, code may be injected via `user_callbacks`,
# which is a `Tuple` of [`GenericCallbacks`](@ref).
ClimateMachine.invoke!(solver_config; user_callbacks = (callback,));

# # Post-processing

# Our solution is stored in the nested dictionary `all_data` whose keys are
# the output interval. The next level keys are the variable names, and the
# values are the values along the grid:

# To get `T` at ``t=0``, we can use `T_at_t_0 = all_data[0]["T"][:]`
@show keys(all_data[0])

# Let's plot the solution:
####

z_scale = 100 # convert from meters to cm
z_key = "z"
z_label = "z [cm]"
#What is get_z - how does this map the grid to z_scale? Just a scalar?
z = get_z(mygrid, z_scale)

output_dir = @__DIR__
export_plot(
    -z,
    all_data,
    ("ν",),
    joinpath(output_dir, "foo.png"),
    z_label,
);
