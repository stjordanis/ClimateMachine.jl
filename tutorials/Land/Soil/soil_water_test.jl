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
using DelimitedFiles
using Parameters

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
include("soil_water_model_test.jl")

# Add water functions
include("Water/water_functions.jl")
include("Water/matric_potential_composable.jl")
include("Water/hydraulic_conductivity_composable.jl")
WF = waterfunctions(
    hydraulic_cond = vanGenuchten{FT}(),
    matric_pot = vanGenuchten{FT}()
)

# Define time variables
const minute = 60
const hour = 60*minute
const day = 24*hour

######
###### 2) Set up system
######
println("2) Set up system...")

# Read in state variables and data. to be specified using a soil parameters struct and treated as an attribute of the model.
mineral_properties = "Clay"
K_sat  = 0.0443 / (3600*100)
porosity = 0.495 # Read in from data base
S_s = 10e-4  # [ m-1] ## Check this value !!

#IC/BC values  - to be specified/calculated via a BC/IC struct, and treated as an attribute of the model.
ν_0 = 0.24
ν_surface = porosity-1e-3
S_l_0 = effective_saturation(porosity, ν_0)
ψ_0 = pressure_head(WF.matric_pot,S_l_0,porosity,S_s,ν_0)
println(ψ_0)
κ_0 = hydraulic_conductivity(WF.hydraulic_cond, K_sat, S_l_0, ψ_0,0.0)


# Load Soil Model in 'm'
m = SoilModelMoisture(
    param_set = my_param_set,
    WF = WF,
    # Define hydraulic conductivity of soil
    initialκ   = (aux) -> κ_0,
    # Define initial soil moisture
    initialν = (state, aux) -> ν_0, #theta_liq_0, # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287,
    surfaceν = (state, aux, t) -> ν_surface, #theta_liq_surface, # [m3/m3] constant flux at surface, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
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
zmax = FT(0);


# Establish a `ClimateMachine` single stack configuration
driver_config = ClimateMachine.SingleStackConfiguration(
    "SoilMoistureModel",
    N_poly,
    nelem_vert,
    zmax,
    my_param_set,
    m;
    zmin = FT(-1),
    numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
);

# Minimum spatial and temporal steps
Δ = min_node_distance(driver_config.grid)
τ = (Δ^2 /K_sat)
dt = 6#0.002*τ #CFL_bound*0.5 # TODO: check if this is a reasonable expression




######
###### 4) Prep ICs, and time-stepper and output configurations
######
println("4) Prep ICs, and time-stepper and output configurations...")
# # Configure a `ClimateMachine` solver.
const timeend = FT(1*day)
const t0 = FT(0)

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

# Define the number of outputs from `t0` to `timeend`
const n_outputs = 5;

# This equates to exports every ceil(Int, timeend/n_outputs) time-step:
const every_x_simulation_time = ceil(Int, timeend / n_outputs);

# Create a nested dictionary to store the solution:
all_data = Dict([k => Dict() for k in 0:n_outputs]...)
# The `ClimateMachine`'s time-steppers provide hooks, or callbacks, which
# allow users to inject code to be executed at specified intervals. In this
# callback, the state and aux variables are collected, combined into a single
# `OrderedDict` and written to a NetCDF file (for each output step `step`).
step = [0];
callback = GenericCallbacks.EveryXSimulationTime(
    every_x_simulation_time,
    solver_config.solver,
) do (init = false)
    t = ODESolvers.gettime(
        solver_config.solver
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
    all_vars = OrderedDict(state_vars..., aux_vars...)
    all_vars["t"]= [t]
    all_data[step[1]] = all_vars

    step[1] += 1
    nothing
end;

# # Solve

# This is the main `ClimateMachine` solver invocation. While users do not have
# access to the time-stepping loop, code may be injected via `user_callbacks`,
# which is a `Tuple` of [`GenericCallbacks`](@ref).
ClimateMachine.invoke!(solver_config; user_callbacks = (callback,));
#Pull out state etc. at final time step
t = ODESolvers.gettime(
    solver_config.solver
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
all_vars = OrderedDict(state_vars..., aux_vars...);
all_vars["t"]= [t]
all_data[n_outputs] = all_vars
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
    z,
    all_data,
    ("ν",),
    joinpath(output_dir, "foo.png"),
    z_label,
);


open("./final_step.txt", "w") do io
    writedlm(io, all_data[n_outputs])
end
