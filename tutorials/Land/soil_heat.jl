# --------------------------------- run CLIMA SOIL MODEL -----------------------
# CLIMA_run_SoilHeat.jl: This model simulates soil heat dynamics for the CliMA model

######
###### 1) Import/Export Needed Functions
######
println("1) Import/Export Needed Functions")

# Load necessary CliMA subroutines
using MPI
using Test
using ClimateMachine
using Logging
using Printf
using NCDatasets
using LinearAlgebra
using OrderedCollections
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.Writers
using ClimateMachine.VTK
using ClimateMachine.Mesh.Elements: interpolationmatrix
using ClimateMachine.DGmethods
using ClimateMachine.DGmethods.NumericalFluxes
using ClimateMachine.MPIStateArrays
using ClimateMachine.GenericCallbacks: EveryXWallTimeSeconds, EveryXSimulationSteps
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers

using Interpolations
using DelimitedFiles

# ENV["GKS_ENCODING"] = "utf-8"

ENV["CLIMA_GPU"] = "false"

# Initialize CliMA
ClimateMachine.init()

FT = Float64

# Change output directory and save plots there
output_dir = joinpath(dirname(dirname(pathof(ClimateMachine))), "output", "land")
mkpath(output_dir)

# Add soil model
include("soil_heat_model.jl")
# Add other functions
include(joinpath("Heat","thermal_properties.jl"))
include(joinpath("Heat","kersten.jl"))
include(joinpath("Heat","heat_capacity.jl"))
include(joinpath("Heat","internal_energy.jl"))
include(joinpath("Heat","temperature_calculator.jl"))

######
###### Include helper and plotting functions (to be refactored/moved into CLIMA src)
######

include("helper_funcs.jl")
include("plotting_funcs.jl")

# Read in real temperature data
Real_Data_vector =  readdlm("tutorials/Land/Heat/T_desert_25N_25E.txt", '\t', FT, '\n')
Real_Data_vector = [Real_Data_vector...]
Real_Data_vector = collect(Real_Data_vector)

# discrete time
# Real_time_data = range(0.00,stop=31536000,length= 8760)
Real_time_data = range(0.00,stop=31536000,length= length(Real_Data_vector))
Real_time_data = collect(Real_time_data)

# return a data structure that is a callable function
Real_continuous_data = TimeContinuousData(Real_time_data, Real_Data_vector)

######
###### 2) Set up domain
######
println("2) Set up domain...")

# Read in state variables and data
theta_liquid_0 = 0.35 # Read in from water model {state.θ}
theta_ice_0 = 0.05 # Read in from water model {state.θi}
mineral_properties = "Sand" # Read in from data base
porosity = 0.5 # Read in from data base

# NOTE: this is using 5 vertical elements, each with a 5th degree polynomial,
# giving an approximate resolution of 5cm
const velems = 0.0:-0.2:-1 # Elements at: [0.0 -0.2 -0.4 -0.6 -0.8 -1.0] (m)
const N = 5 # Order of polynomial function between each element

# Set domain using Stacked Brick
grid = SingleStackGrid(MPI, velems, N, FT, Array)

# Define SoilModel struct
m = SoilModel(
    # Define heat capacity of soil
     ρc = (state, aux, t) -> heat_capacity(mineral_properties,porosity,theta_liquid_0,theta_ice_0 ), # state.θ,state.θi
    # ρc = (state, aux, t) ->  aux.z > -0.5 ? 2.49e6 : 2.61e6,

    # Define thermal conductivity of soil
     κ   = (state, aux, t) ->   thermal_properties(mineral_properties,theta_liquid_0,theta_ice_0), # state.θ,state.θi
    # κ  = (state, aux, t) ->  aux.z > -0.5 ? 2.42 : 1.17,

    # Define initial temperature of soil
    initialT = (aux, t) -> (273.15 + 12.0),

    # Define surface boundary condition
    #surfaceT = (state, aux, t) -> (273.15 + 12.0) + 0.5*10.0 * sinpi(2*(t/(60*60)-8)/24)
    #surfaceT = (state, aux, t) -> (273.15 + 12.0) + 250*(1/sqrt(2*pi*10^2))*exp( -((t/(60*60)-24)^2)/(2*10^2) )
    surfaceT = (state, aux, t) -> Real_continuous_data(t) # replace with T_data
)

# Set up DG scheme
dg = DGModel( #
  m, # "PDE part"
  grid,
  CentralNumericalFluxFirstOrder(), # penalty terms for discretizations
  CentralNumericalFluxSecondOrder(),
  CentralNumericalFluxGradient())

# Minimum spatial and temporal steps
Δ = min_node_distance(grid)
CFL_bound = (Δ^2 / (2 * 2.42/2.49e6))
dt = CFL_bound*0.5 # TODO: provide a "default" timestep based on  Δx,Δy,Δz


######
###### 3) Define variables for simulation
######
println("3) Define variables for simulation...")

# Define time variables
const minute = 60
const hour = 60*minute
const day = 24*hour
# const timeend = 1*minute
# const n_outputs = 25
const timeend = 5*hour

# Output frequency:
# const every_x_simulation_time = ceil(Int, timeend/n_outputs)
const every_x_simulation_time = 1*hour


######
###### 4) Prep ICs, and time-stepper and output configurations
######
println("4) Prep ICs, and time-stepper and output configurations...")

# state variable
Q = init_ode_state(dg, Float64(0))

# initialize ODE solver
lsrk = LSRK54CarpenterKennedy(dg, Q; dt = dt, t0 = 0)

mkpath(output_dir)

dims = OrderedDict("z" => collect(get_z(grid, 100)))
# run for 8 days (hours?) to get to steady state

output_data = DataFile(joinpath(output_dir, "output_data"))

step = [0]
stcb = GenericCallbacks.EveryXSimulationTime(every_x_simulation_time, lsrk) do (init = false)
  state_vars = get_vars_from_stack(grid, Q, m, vars_state_conservative) # ; exclude=["θi"])
  aux_vars = get_vars_from_stack(grid, dg.state_auxiliary, m, vars_state_auxiliary; exclude=["z"])
  all_vars = OrderedDict(state_vars..., aux_vars...)
  write_data(NetCDFWriter(), output_data(step[1]), dims, all_vars, gettime(lsrk))
  step[1]+=1
  nothing
end


######
###### 5) Solve the equations
######
println("5) Solve the equations...")

solve!(Q, lsrk; timeend=timeend, callbacks=(stcb,))

#####
##### 6) Post-processing
#####
println("6) Post-processing...")

all_data = collect_data(output_data, step[1])

# To get "T" at timestep 0:
# all_data[0]["T"][:]
