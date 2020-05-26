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
const timeend = hour

# Output frequency:
# const every_x_simulation_time = ceil(Int, timeend/n_outputs)
const every_x_simulation_time = 10*minute

######
###### 3) # Add soil model and other functions
######
println("3) Add soil model and other functions...")


function heaviside(t) # Should I be using an approximation of a heaviside function so that it is continuous
   0.5 * (sign(t) + 1)
end
t1=1*hour
t2=t1+10*day

include("soil_heat_model_unit_test.jl")
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

#### Energy conservation check

#using Pkg
#Pkg.add("SymPy")
using SymPy
t=0.0:-0.1:-1

E=rand(Float64, length(all_data)-1)
for i=1:length(all_data)
   E[i]=all_data[1]["int.a"][end]
end

#E=rand(Float64, length(all_data)-1)
#for i=2:length(all_data)-1
#  E[i-1]=all_data[i]["int.a"][end]-all_data[1]["int.a"][end]
#end

#δE_analytical=rand(Float64, length(all_data)-2)
##δE_analytical_remain=rand(Float64, length(all_data)-2)
#for i=1:length(all_data)-2
#  δE_analytical[i]=10*2.42*every_x_simulation_time*i
#  #δE_analytical_remain[i]=10*2.42*every_x_simulation_time*i/dt-round(10*2.42*every_x_simulation_time*i/dt)
#end

δE_analytical=rand(Float64, length(all_data))
#δE_analytical_remain=rand(Float64, length(all_data)-2)
for i=1:length(all_data)
   δE_analytical[i]=10*2.42*every_x_simulation_time*(i)+2.49e6*273.15
   #δE_analytical_remain[i]=10*2.42*every_x_simulation_time*i/dt-round(10*2.42*every_x_simulation_time*i/dt)
end

## error over time
#error=E-δE_analytical
#rel_error=error[end]/(E[end]+all_data[1]["int.a"][end])
#rel_error2=error[end]/(δE_analytical[end])
#using Plots
#
x = range(0, stop=timeend, length=length(E))
display(plot(x,[E δE_analytical],
xlabel = "Time [s]",
ylabel = "Energy [J]",
label = ["E" "δE_analytical"]))

# L2 Norm, report
# run with smaller time steps, make plot of error as function of dt
#||v||2 = sqrt(a1^2 + a2^2 + a3^2)

##### Analytical vs numerical solution
#using Pkg
#Pkg.add("SymPy")
#using SymPy
#x = Sym("x")
#k = Sym("k")
#t = Sym("t")
#
#minute=60
#hour=60*minute
#day=24*hour
#t1=1*hour
#t2=t1+10*day
#
#alpha=2.42/2.49e6
#L=-1
#A=10
#Ap=0
#B=0
#Bp=0
#p=0
#lambda=k*pi/L
#f=273.15
#
#g=f-A*x-(B-A)*x^2/(2*L)
#q=p-Ap*x-(Bp-Ap)*x^2/(2*L)+alpha/L*(B-A)
#a0=(1/L)*-16589/60 #int(g,x,0,L)
#an=(2/L)*-(10*sinpi(k) - 10*pi*k + (5563*k^2*pi^2*sinpi(k))/20)/(k^3*pi^3) #int(g*cos(lambda*x),x,0,L)
#A0=(1/L)*-1434252872558301/147573952589676412928 #int(q,x,0,L)
#An=(2/L)*-(1434252872558301*sinpi(k))/(147573952589676412928*k*pi)#int(q*cos(lambda*x),x,0,L)
#u1p=A.*x+(B-A).*x^2/(2*L)
#s1=an.*exp(-alpha*t*lambda^2)*cos(lambda*x)
#u3p=-(1434252872558301*t)/147573952589676412928 #int(A0,s,0,t)
#s2=cos(lambda*x)*-(91792183843731264*sinpi(k)*(exp(-(4589609192186563*k^2*t*pi^2)/4722366482869645213696) - 1))/(4589609192186563*k^3*pi^3) #int(An*exp(-alpha*lambda^2*(t-s)),s,0,t)
#
#using Plots
#
##for t=hour #0:10*minute:hour
#t=hour
#x=0:L:1000
#u1=eval(u1p,x)
#u2=sum(k-> s1, 1:1000)  #u2=eval(symsum(s1,k,1,1000))
#u3=eval(u3p)
#u4=sum(k-> s2, 1:1000)  #u4=eval(symsum(s2,k,1,1000))
#u=u1+a0+u2+u3+u4 #,'Color',c(i,:),'LineWidth',3)
##end
#
#display(plot(x,u,
#  xlabel = "Time [s]",
#ylabel = "Energy [J]",
#label = ["E" "δE_analytical"]))
