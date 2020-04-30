# --------------------------------- run CLIMA SOIL MODEL -----------------------
# CLIMA_run_SoilHeat.jl: This model simulates soil heat dynamics for the CliMA model

######
###### 1) Import/Export Needed Functions
######
println("1) Import/Export Needed Functions")

# Load necessary CliMA subroutines
using MPI
using Test
using CLIMA
using Logging
using Printf
using NCDatasets
using LinearAlgebra
using OrderedCollections
using CLIMA.Mesh.Topologies
using CLIMA.Mesh.Grids
using CLIMA.Writers
using CLIMA.VTK
using CLIMA.Mesh.Elements: interpolationmatrix
using CLIMA.DGmethods
using CLIMA.DGmethods.NumericalFluxes
using CLIMA.MPIStateArrays
using CLIMA.GenericCallbacks: EveryXWallTimeSeconds, EveryXSimulationSteps
using CLIMA.GenericCallbacks
using CLIMA.ODESolvers

using Interpolations
using DelimitedFiles

# ENV["GKS_ENCODING"] = "utf-8"

ENV["CLIMA_GPU"] = "false"

# Initialize CliMA
CLIMA.init()

FT = Float64

# Change output directory and save plots there
output_dir = joinpath(dirname(dirname(pathof(CLIMA))), "output", "land")
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
const timeend = 18*day

# Output frequency:
# const every_x_simulation_time = ceil(Int, timeend/n_outputs)
const every_x_simulation_time = 1*day

######
###### 3) # Add soil model and other functions
######
println("3) Add soil model and other functions...")

include("CLIMA_SoilHeat_unit_test.jl")
#include("thermal_properties.jl")
#include("kersten.jl")

######
###### Include helper and plotting functions (to be refactored/moved into CLIMA src)
######

function get_z(grid)
    # TODO: this currently uses some internals: provide a better way to do this
    return reshape(grid.vgeo[(1:(N+1)^2:(N+1)^3),CLIMA.Mesh.Grids.vgeoid.x3id,:],:)*100
end
include(joinpath("..","helper_funcs.jl"))
include(joinpath("..","plotting_funcs.jl"))

# Read in real temperature data
Real_Data_vector =  readdlm("examples/Land/Heat/T_desert_25N_25E.txt", '\t', FT, '\n')
Real_Data_vector = [Real_Data_vector...]
Real_Data_vector = collect(Real_Data_vector)

# discrete time
# Real_time_data = range(0.00,stop=31536000,length= 8760)
Real_time_data = range(0.00,stop=31536000,length= length(Real_Data_vector))
Real_time_data = collect(Real_time_data)

# return a data structure that is a callable function
Real_continuous_data = TimeContinuousData(Real_time_data, Real_Data_vector)

######
###### 4) Set up domain
######
println("4) Set up domain...")

# NOTE: this is using 5 vertical elements, each with a 5th degree polynomial
# giving an approximate resolution of 5cm (true when elements at: [0.0 -0.2 -0.4 -0.6 -0.8 -1.0] and 5th degree polynomial)
const velems = 0.0:-0.1:-1 # Elements at: [0.0 -0.2 -0.4 -0.6 -0.8 -1.0] (m)
const N = 4 # Order of polynomial function between each element

# Set domain using Stached Brick Topology
topl = StackedBrickTopology(MPI.COMM_WORLD, (0.0:1,0.0:1,velems);
    periodicity = (true,true,false),
    boundary=((0,0),(0,0),(1,2)))

# Set up grid
grid = DiscontinuousSpectralElementGrid(topl, FloatType = FT, DeviceArray = Array, polynomialorder = N)

# Define thermal conductivity
#κ_sand = (thermal_properties("Sand",0.35,0.05 ))
#κ_clay = (thermal_properties("Clay",0.35,0.05 ))
#κ_other = (thermal_properties("Other",0.35,0.05 ))

function heaviside(t) # Should I be using an approximation of a heaviside function so that it is continuous
   0.5 * (sign(t) + 1)
end
t1=1*hour
#t2=t1+10*day

# Define SoilModel struct
m = SoilModel(
    ρc = (state, aux, t) -> 2.49e6, # aux.z > -0.5 ? 2.49e6 : 2.61e6,
    κ  = (state, aux, t) -> 2.42, # aux.z > -0.5 ? κ_sand : κ_clay,
    initialT = (aux) -> 273.15, #  (273.15 + 2 + 5*exp(-(aux.z-0.5)^2/(2*(0.2)^2)))),
    surfaceT = (state, aux, t) -> 273.15+10*heaviside(t-t1) #-10*heaviside(t-t2) #(273.15 + 15.0 + 0.5*10.0 * sinpi(2*(t/(60*60)-8)/24)) # replace with T_data
    #surfaceT = (state, aux, t) -> (273.15 + 12.0) + 250*(1/sqrt(2*pi*10^2))*exp( -((t/(60*60)-24)^2)/(2*10^2) ) # replace with T_data
    #surfaceT = (state, aux, t) -> Real_continuous_data(t) # replace with T_data
)

# Set up DG scheme
dg = DGModel( #
  m, # "PDE part"
  grid,
  CentralNumericalFluxNonDiffusive(), # penalty terms for discretizations
  CentralNumericalFluxDiffusive(),
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

# Plot initial state
p = get_plot(grid, Q, dg.auxstate, 0)
export_plots(p, joinpath(output_dir, "initial_state.png"))

mkpath(output_dir)

plots = []
dims = OrderedDict("z" => collect(get_z(grid)))
# run for 8 days (hours?) to get to steady state

output_data = DataFile(joinpath(output_dir, "output_data"))

step = [0]
stcb = GenericCallbacks.EveryXSimulationTime(every_x_simulation_time, lsrk) do (init = false)
  state_vars = get_vars_from_stack(grid, Q, m, vars_state)
  aux_vars = get_vars_from_stack(grid, dg.auxstate, m, vars_aux; exclude=["z"])
  integral_vars = get_vars_from_stack(grid, Q, m, vars_integral)
  all_vars = OrderedDict(state_vars..., aux_vars..., integral_vars...)
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

# OLD:
# plots = []
# for num in keys(all_data)
#   Tg = all_data[num]["Tg"]
#   p = plot(Tg, gridg, ylabel="depth (cm) at t=$(t)", xlabel="T (°K)", yticks=-100:20:0, xlimits=(263.15,303.15), legend=false)
#   push!(plots, plot)

# export_plots(plots, joinpath(output_dir, "state_over_time.png"))
# a function for performing interpolation on the DG grid
# TODO: use CLIMA interpolation once available

# t_plot = 24*7 # How many time steps to plot?
# t_plot = 1 # How many time steps to plot?
# Zprofile = -0.995:0.01:0 # needs to be in sorted order for contour function
# Tprofile = zeros(length(Zprofile),t_plot)
# hours = 0.5:1:t_plot

# solve! should occur only once, not in a loop
# for (i,h) in enumerate(hours)
#    t = solve!(Q, lsrk; timeend=0*day+h*hour)
#    Tprofile[:,i] = (interpolate_grid(grid, dg.auxstate, CLIMA.Mesh.Grids.vgeoid.x3id, Zprofile))
# end

# plot_contour(hours, Zprofile, Tprofile, t_plot, filename) = nothing


# plot_contour(hours, Zprofile, Tprofile, t_plot, joinpath(output_dir, "contour_T.png"))
