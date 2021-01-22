using Test
using MPI

using ClimateMachine

# To test coupling
using ClimateMachine.Coupling

# To create meshes (borrowed from Ocean for now!)
using ClimateMachine.Ocean.Domains

# To setup some callbacks
using ClimateMachine.GenericCallbacks

# To invoke timestepper
using ClimateMachine.ODESolvers

ClimateMachine.init()

# Use toy balance law for now
include("CplTestingBL.jl")

# Make some meshes covering same space laterally.
domainA = RectangularDomain(
    Ne = (10, 10, 5),
    Np = 4,
    x = (0, 1e6),
    y = (0, 1e6),
    z = (0, 1e5),
    periodicity = (true, true, true),
)
domainO = RectangularDomain(
    Ne = (10, 10, 4),
    Np = 4,
    x = (0, 1e6),
    y = (0, 1e6),
    z = (-4e3, 0),
    periodicity = (true, true, true),
)
domainL = RectangularDomain(
    Ne = (10, 10, 1),
    Np = 4,
    x = (0, 1e6),
    y = (0, 1e6),
    z = (0, 1),
    periodicity = (true, true, true),
)

# Create 3 components - one on each domain, for now all are instances
# of the same balance law
mA=Coupling.CplTestModel(;domain=domainA,BL_module=CplTestingBL, nsteps=5)
mO=Coupling.CplTestModel(;domain=domainO,BL_module=CplTestingBL, nsteps=2)
#mL=Coupling.CplTestModel(;domain=domainL,BL_module=CplTestingBL)

function postatmos(_)
    println("Atmos export fill callback")
    mO.discretization.state_auxiliary.boundary_in[...] .= mA.discretization.state_auxiliary.boundary_out[...]

end
function postocean(_)
    println("Ocean export fill callback")
    mA.discretization.state_auxiliary.boundary_in[...] .= mO.discretization.state_auxiliary.boundary_out[...]
end

# Instantiate a coupled timestepper that steps forward the components and
# implements mapings between components export bondary states and
# other components imports.
component_list=(atmosphere=mA,ocean=mO,)
callback_list =(atmosphere=postatmos,ocean=postocean,)
cC=Coupling.CplSolver(component_list=component_list,
                      callback_list=callback_list,
                      coupling_dt=5.,t0=0.)

# Invoke solve! with coupled timestepper and callback list.
solve!(nothing,cC;numberofsteps=2)
