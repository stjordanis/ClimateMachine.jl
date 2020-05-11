# --------------------------------- run CLIMA SOIL MODEL -----------------------
# CLIMA_run_SoilWater.jl: This model simulates soil water dynamics for the CliMA model

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

# Add soil moisture model
include("CLIMA_SoilWater.jl")
# Add other functions
include("soil_water_properties.jl")
include("frozen_impedence_factor.jl")
include("temperature_dependence.jl")
include("matric_potential.jl")
include("pressure_head.jl")
include("hydraulic_head.jl")
include("effective_saturation.jl")
include("augmented_liquid.jl")
include("calculate_frozen_water.jl")


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
mineral_properties = "Sand"
soil_T = 265 # Read in from heat model {aux.T}
soil_Tref = 282.42 # Soil reference temperature: annual mean temperature of site
theta_liq_0 = 0.25 # Read in from water model {state.θ}
theta_liq_surface = 0.25 # Read in from water model {state.θ}
theta_ice_0 = 0.02 # Read in from water model {state.θi}
h_0 = -3 # Read in from water model {state.θ}
ψ_0 = -1 # Soil pressure head {aux.h}
porosity = 0.8 # Read in from data base
S_s = 10e-4  # [ m-1]
flag = "van Genuchten" # "van Genuchten" , "Brooks and Corey"


# NOTE: this is using 5 vertical elements, each with a 5th degree polynomial,
# giving an approximate resolution of 5cm
const velems = 0.0:-0.2:-1 # Elements at: [0.0 -0.2 -0.4 -0.6 -0.8 -1.0] (m)
const N = 5 # Order of polynomial function between each element

# Set domain using Stached Brick Topology
topl = StackedBrickTopology(MPI.COMM_WORLD, (0.0:1,0.0:1,velems);
    periodicity = (true,true,false),
    boundary=((0,0),(0,0),(1,2)))
grid = DiscontinuousSpectralElementGrid(topl, FloatType = Float64, DeviceArray = Array, polynomialorder = N)

# Load Soil Model in 'm'
m = SoilModelMoisture(
     # Define hydraulic conductivity of soil
     K_s   = (state, aux, t) ->   soil_water_properties(mineral_properties,soil_T,soil_Tref,state.θ,state.θi,porosity,aux.ψ,S_s,flag), #aux.T,state.θ,state.θi,aux.h 
    # K_s  = (state, aux, t) -> (1e-3*(0.34/(60*60))*1.175e6/((1.175e6+abs.(aux.h-aux.z)^4.74))), #(0.34)
    
    # Define initial soil moisture
    initialθ = (aux, t) -> theta_liq_0, # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287,
    surfaceθ = (state, aux, t) -> theta_liq_surface, # [m3/m3] constant flux at surface, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
    
    # Define initial and boundary condition parameters
    initialh = (aux, t) -> h_0, #100- aux.z  # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
    
    # Define initial and boundary condition parameters
    initialψ = (aux, t) -> h_0 - aux.z, # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
    
    # Define initial and boundary condition parameters
    initialθi = (aux, t) -> theta_ice_0  #267 # [m3/m3] constant flux at surface, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287

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
###### 3) Define variables for simulation
######
println("3) Define variables for simulation...")

# Define time variables
const minute = 60
const hour = 60*minute
const day = 24*hour
# const timeend = 1*minute
# const n_outputs = 25
const timeend = 1*hour

# Output frequency:
# const every_x_simulation_time = ceil(Int, timeend/n_outputs)
const every_x_simulation_time = 10*minute


######
###### 4) Prep ICs, and time-stepper and output configurations
######
println("4) Prep ICs, and time-stepper and output configurations...")


# state variable
Q = init_ode_state(dg, Float64(0))

# initialize ODE solver
lsrk = LSRK54CarpenterKennedy(dg, Q; dt = dt, t0 = 0)

# Plot initial state
p = get_plot(grid, Q, dg.auxstate, 0)
export_plots(p, joinpath(output_dir, "initial_state_Water.png"))

mkpath(output_dir)

plots = []
dims = OrderedDict("z" => collect(get_z(grid, 100)))
# run for 8 days (hours?) to get to steady state

output_data = DataFile(joinpath(output_dir, "output_data_Water"))

step = [0]
stcb = GenericCallbacks.EveryXSimulationTime(every_x_simulation_time, lsrk) do (init = false)
  state_vars = get_vars_from_stack(grid, Q, m, vars_state) #; exclude=["θi"])
  aux_vars = get_vars_from_stack(grid, dg.auxstate, m, vars_aux; exclude=["z"])
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



# OLD ------------------------ 5) Run model for many time steps and plot ---------------------------------------

# a function for performing interpolation on the DG grid
# TODO: use CLIMA interpolation once available
#function interpolate(grid, auxstate, Zprofile)
#    P = zeros(size(Zprofile))
#    nelems = size(grid.vgeo, 3)
#    for elem in 1:nelems
#        G = grid.vgeo[(1:(N+1)^2:(N+1)^3),CLIMA.Mesh.Grids.vgeoid.x3id,elem]
#        I = minimum(G) .< Zprofile .<= maximum(G)
#        M = interpolationmatrix(G, Zprofile[I])
#        P[I] .= M*auxstate.data[(1:(N+1)^2:(N+1)^3),2,elem]
#    end
#    return P
#end

#t_plot = 24*4 # How many time steps to plot?
#Zprofile = -0.995:0.01:0 # needs to be in sorted order for contour function
#Tprofile = zeros(length(Zprofile),t_plot)
#hours = 0.5:1:t_plot

#for (i,h) in enumerate(hours)
#    t = solve!(Q, lsrk; timeend=day+h*hour)
#    Tprofile[:,i] = (interpolate(grid, dg.auxstate, Zprofile))
#end

#contour(hours, Zprofile.*100, Tprofile,
#    levels=0:1, xticks=0:4:t_plot, xlimits=(0,t_plot),
#    xlabel="Time of day (hours)", ylabel="Soil depth (cm)", title="Volumetric water content (m3/m3)")

#savefig(joinpath(output_dir, "contour.png"))
