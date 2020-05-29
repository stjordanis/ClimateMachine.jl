# --------------------------------- run CLIMA SOIL MODEL -----------------------
# CLIMA_run_SoilWater.jl: This model simulates soil water dynamics for the CliMA model

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

# Add soil moisture model
include("soil_water_model.jl")

# Add water functions
#include("Water/soil_water_properties.jl")
#include("Water/frozen_impedence_factor.jl")
#include("Water/temperature_dependence.jl")
include("Water/matric_potential.jl")
include("Water/pressure_head.jl")
include("Water/hydraulic_head.jl")
include("Water/effective_saturation.jl")
#include("Water/augmented_liquid.jl")
#include("Water/calculate_frozen_water.jl")
#include("Water/heaviside.jl")

######
###### Include helper and plotting functions (to be refactored/moved into CLIMA src)
######

include(joinpath("..","helper_funcs.jl"))
include(joinpath("..","plotting_funcs.jl"))

# # Read in real temperature data
# Real_Data_vector =  readdlm("examples/Land/Heat/T_desert_25N_25E.txt", '\t', FT, '\n')
# Real_Data_vector = [Real_Data_vector...]
# Real_Data_vector = collect(Real_Data_vector)

# # discrete time
# # Real_time_data = range(0.00,stop=31536000,length= 8760)
# Real_time_data = range(0.00,stop=31536000,length= length(Real_Data_vector))
# Real_time_data = collect(Real_time_data)

# # return a data structure that is a callable function
# Real_continuous_data = TimeContinuousData(Real_time_data, Real_Data_vector)

######
###### 2) Set up domain
######
println("2) Set up domain...")

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
# NOTE: this is using 5 vertical elements, each with a 5th degree polynomial,
# giving an approximate resolution of 5cm
const velems = 0.0:-0.1:-1 # Elements at: [0.0 -0.2 -0.4 -0.6 -0.8 -1.0] (m)
const N = 4 # Order of polynomial function between each element

# Set domain using Stacked Brick
grid = SingleStackGrid(MPI, velems, N, FT, Array)

# Load Soil Model in 'm'
m = SoilModelMoisture(
     # Define hydraulic conductivity of soil
    K_s   = (state, aux, t) ->  K_sat*124.6/(124.6 + abs(aux.h-aux.z)^1.77), #soil_water_properties(mineral_properties,soil_T,soil_Tref,state.θ,state.θi,porosity,aux.ψ,S_s,flag), #aux.T,state.θ,state.θi,aux.h
    # Define initial soil moisture
    initialν = (state, aux) -> ν_0, #theta_liq_0, # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287,
    surfaceν = (state, aux, t) -> ν_surface, #theta_liq_surface, # [m3/m3] constant flux at surface, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
    # Define initial and boundary condition parameters
    #initialS_l = (aux) -> S_l_0,
    #initialψ_m = (aux) -> ψ_m_0,
    #initialψ = (aux) -> ψ_0,
    initialh = (aux) -> aux.z + ψ_0 # [m3/m3] constant water content in soil
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
dt = 20 #CFL_bound*0.5 # TODO: provide a "default" timestep based on  Δx,Δy,Δz

######
###### 3) Define variables for simulation
######
println("3) Define variables for simulation...") # move up

# Define time variables
const minute = 60
const hour = 60*minute
const day = 24*hour
# const timeend = 1*minute
# const n_outputs = 25
const timeend = 24*day

# Output frequency:
# const every_x_simulation_time = ceil(Int, timeend/n_outputs)
const every_x_simulation_time = 1*day


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
  state_vars = get_vars_from_stack(grid, Q, m, vars_state_conservative)
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

#using Plots
#
#display(plot(x,u,
# xlabel = "Time [s]",
#ylabel = "Energy [J]",
#label = ["E" "δE_analytical"]))