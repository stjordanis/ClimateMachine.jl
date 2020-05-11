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


######
###### 2) Define variables for simulation
######
println("2) Define variables for simulation...")

# Define time variables
const minute = 60
const hour = 60*minute
const day = 24*hour
# const timeend = 1*minute
# const n_outputs = 25
const timeend = 5*day

# Output frequency:
# const every_x_simulation_time = ceil(Int, timeend/n_outputs)
const every_x_simulation_time = 6*hour

######
###### 3) # Add soil model and other functions
######
println("3) Add soil model and other functions...")


function heaviside(t) # Should I be using an approximation of a heaviside function so that it is continuous
   0.5 * (sign(t) + 1)
end
t1=1*hour
t2=t1+10*day

include("CLIMA_SoilHeat_unit_test.jl")
#include("thermal_properties.jl")
#include("kersten.jl")

######
###### Include helper and plotting functions (to be refactored/moved into CLIMA src)
######

include(joinpath("..","helper_funcs.jl"))
include(joinpath("..","plotting_funcs.jl"))

# Read in real temperature data
#Real_Data_vector =  readdlm("examples/Land/Heat/T_desert_25N_25E.txt", '\t', FT, '\n')
#Real_Data_vector = [Real_Data_vector...]
#Real_Data_vector = collect(Real_Data_vector)

# discrete time
# Real_time_data = range(0.00,stop=31536000,length= 8760)
#Real_time_data = range(0.00,stop=31536000,length= length(Real_Data_vector))
#Real_time_data = collect(Real_time_data)

# return a data structure that is a callable function
#Real_continuous_data = TimeContinuousData(Real_time_data, Real_Data_vector)

######
###### 4) Set up domain
######
println("4) Set up domain...")

# NOTE: this is using 5 vertical elements, each with a 5th degree polynomial
# giving an approximate resolution of 5cm (true when elements at: [0.0 -0.2 -0.4 -0.6 -0.8 -1.0] and 5th degree polynomial)
const velems = 0.0:-0.1:-1 # Elements at: [0.0 -0.2 -0.4 -0.6 -0.8 -1.0] (m)
const N = 4 # Order of polynomial function between each element

# Set domain using Stacked Brick
grid = SingleStackGrid(MPI, velems, N, FT, Array)

# Define thermal conductivity
#κ_sand = (thermal_properties("Sand",0.35,0.05 ))
#κ_clay = (thermal_properties("Clay",0.35,0.05 ))
#κ_other = (thermal_properties("Other",0.35,0.05 ))

# Define SoilModel struct
m = SoilModel(
    ρc = (state, aux, t) -> 2.49e6, # aux.z > -0.5 ? 2.49e6 : 2.61e6,
    κ  = (state, aux, t) -> 2.42, # aux.z > -0.5 ? κ_sand : κ_clay,
    initialT = (aux) -> 273.15, # aux.z - (aux.z)^2,  #  (273.15 + 2 + 5*exp(-(aux.z-0.5)^2/(2*(0.2)^2)))),
    surfaceT = (state, aux, t) -> 273.15 #+10*heaviside(t-t1) #-10*heaviside(t-t2) #(273.15 + 15.0 + 0.5*10.0 * sinpi(2*(t/(60*60)-8)/24)) # replace with T_data
    #surfaceT = (state, aux, t) -> (273.15 + 12.0) + 250*(1/sqrt(2*pi*10^2))*exp( -((t/(60*60)-24)^2)/(2*10^2) ) # replace with T_data
    #surfaceT = (state, aux, t) -> Real_continuous_data(t) # replace with T_data
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
###### 5) Prep ICs, and time-stepper and output configurations
######
println("5) Prep ICs, and time-stepper and output configurations...")

# state variable
Q = init_ode_state(dg, FT(0))

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
println("6) Solve the equations...")

solve!(Q, lsrk; timeend=timeend, callbacks=(stcb,))

#####
##### 6) Post-processing
#####
println("7) Post-processing...")

all_data = collect_data(output_data, step[1])

#####
##### 7) Check if analytical solution matches numerical solution
#####
#alpha=2.42/2.49e6; %Change
#L=-1;
#lambda=(2*k-1)*pi/(2*L);
#t0=3600;
#u0=275+10*heaviside(t-t0);
#g=275;
#bk1=2/L*int((g-275)*sin(lambda*x),x,0,L);
#bk2=-10*heaviside(t-t0)*2/L*int(sin(lambda*x),x,0,L);
#bk=bk1+bk2;
#f=bk.*exp(-(alpha)*(lambda^2)*t)*sin(lambda*x);
#figure (1)
#for t=1:3600*24:3600*24*18
#    x=linspace(0,L,1000);
#    u=u0+eval(symsum(f,k,1,1000));
#    u=eval(u);
#    plot(u,x)
#    hold on
#end
#axis([270 290 -1 0])
#ylabel('Depth [m]')
#xlabel('Temperature [K]'


## Check if Gaussian bump is diffusing at right speed

# To get "T" at timestep 0:
# all_data[0]["T"][:]

