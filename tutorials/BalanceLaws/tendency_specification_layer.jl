# # A functional tendency specification layer

# In the balance law (mutating) functions, where we specify fluxes and sources,

# - [`flux_first_order!`](@ref ClimateMachine.BalanceLaws.flux_first_order!)
# - [`flux_second_order!`](@ref ClimateMachine.BalanceLaws.flux_second_order!)
# - and [`source!`](@ref ClimateMachine.BalanceLaws.source!),

# an additional (functional) tendency specification
# layer can be placed on-top that has several nice
# properties. The functional layer:

# - Separates tendency definitions from which tendencies are included in a particular model.
# - Reduces duplicate implementations of tendency definitions (e.g., in optional submodel variants)
# - Allows a more flexible combination of tendencies
# - Allows a simple way to loop over all tendencies for all prognostic variables and recover
#   _each_ flux / source term. This will allow us a simple way to evaluate, for example, the energy budget.

# ## Used modules / imports

# Make running locally easier from ClimateMachine.jl/:
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end

# First, using necessary modules:
using ClimateMachine.BalanceLaws
using ClimateMachine.VariableTemplates
using StaticArrays, Test

# Import methods to overload
import ClimateMachine.BalanceLaws: prognostic_vars, eq_tends, flux

# ## Define a balance law

# Here, we define a simple balance law:

struct MyBalanceLaw <: BalanceLaw end

# ## Define prognostic variable types

# Here, we'll define some prognostic variable types,
# by sub-typing `PrognosticVariable`, for mass and energy:
struct Mass <: PrognosticVariable end
struct Energy <: PrognosticVariable end

# Define [`prognostic_vars`](@ref ClimateMachine.BalanceLaws.prognostic_vars),
# which returns _all_ prognostic variables
prognostic_vars(::MyBalanceLaw) = (Mass(), Energy());

# ## Define some tendencies

# Tendency definitions can be shared across all
# prognostic variables:
struct Advection{PV} <: TendencyDef{Flux{FirstOrder}, PV} end
struct Source1{PV} <: TendencyDef{Source, PV} end
struct Source2{PV} <: TendencyDef{Source, PV} end

# We can also permit tendency definitions to be defined for
# only the relevant prognostic variables: here, we restrict
# the second order flux, Diffusion, to be defined only for Energy.
# Doing this means that if we try to create `Diffusion{Mass}()`,
# an error will be thrown.
struct Diffusion{PV <: Energy} <: TendencyDef{Flux{SecondOrder}, PV} end

# Define [`eq_tends`](@ref ClimateMachine.BalanceLaws.eq_tends),
# which returns a tuple of tendency definitions (those sub-typed
# by [`TendencyDef`](@ref ClimateMachine.BalanceLaws.TendencyDef)),
# given
#  - the prognostic variable
#  - the model (balance law)
#  - the tendency type ([`Flux`](@ref ClimateMachine.BalanceLaws.Flux) or
#    [`Source`](@ref ClimateMachine.BalanceLaws.Source))
#! format: off
eq_tends(pv::PV, ::MyBalanceLaw, ::Flux{FirstOrder}) where {PV <: Mass} = (Advection{PV}(),);
eq_tends(pv::PV, ::MyBalanceLaw, ::Flux{FirstOrder}) where {PV <: Energy} = (Advection{PV}(),);
eq_tends(pv::PV, ::MyBalanceLaw, ::Flux{SecondOrder}) where {PV <: Mass} = ();
eq_tends(pv::PV, ::MyBalanceLaw, ::Flux{SecondOrder}) where {PV <: Energy} = (Diffusion{PV}(),);
eq_tends(pv::PV, ::MyBalanceLaw, ::Source) where {PV <: Mass} = (Source1{PV}(), Source2{PV}());
eq_tends(pv::PV, ::MyBalanceLaw, ::Source) where {PV <: Energy} = (Source1{PV}(), Source2{PV}());
#! format: on

# ## Testing `prognostic_vars` `eq_tends`

# To test that `prognostic_vars` and `eq_tends` were
# implemented correctly, we'll create a balance law
# instance and call [`show_tendencies`](@ref ClimateMachine.BalanceLaws.show_tendencies),
# to make sure that the tendency table is accurate.

bl = MyBalanceLaw()
show_tendencies(bl; table_complete = true)

# The table looks correct. Now we're ready to
# add the specification layer.

# ## Adding the tendency specification layer

# For the purpose of this tutorial, we'll only focus
# on adding the layer to the first order flux, since
# doing so for the second order flux and source
# functions follow the same exact pattern. In other words,
# we'll add a layer that tests the `Flux{FirstOrder}` column
# in the table above. First, we'll define individual
# [`flux`](@ref ClimateMachine.BalanceLaws.flux) kernels:
flux(::Advection{Mass}, bl::MyBalanceLaw, args) = args.state.ρ;
flux(::Advection{Energy}, bl::MyBalanceLaw, args) = args.state.ρe;

# Define `flux_first_order!` and utilize `eq_tends`
function flux_first_order!(
    bl::MyBalanceLaw,
    flx::Grad,
    state::Vars,
    aux,
    t,
    direction,
)

    vec_pad = SVector(1, 1, 1)
    tend = Flux{FirstOrder}()
    args = (; state, aux, t, direction)

    ## `Σfluxes(eq_tends(Mass(), bl, tend), bl, args)` calls
    ## `flux(::Advection{Mass}, ...)` defined above:
    flx.ρ = Σfluxes(eq_tends(Mass(), bl, tend), bl, args) .* vec_pad

    ## `Σfluxes(eq_tends(Energy(), bl, tend), bl, args)` calls
    ## `flux(::Advection{Energy}, ...)` defined above:
    flx.ρe = Σfluxes(eq_tends(Energy(), bl, tend), bl, args) .* vec_pad
    return nothing
end;

# ## Testing the tendency specification layer

# Now, let's test `flux_first_order!` we need to initialize
# some dummy data to call it first:

FT = Float64; # float type
aux = (); # auxiliary fields
t = 0.0; # time
direction = nothing; # Direction

state = Vars{@vars(ρ::FT, ρe::FT)}([1, 2]);
flx = Grad{@vars(ρ::FT, ρe::FT)}(zeros(MArray{Tuple{3, 2}, FT}));

# call `flux_first_order!`
flux_first_order!(bl, flx, state, aux, t, direction);

# Test that `flx` has been properly mutated:
@testset "Test results" begin
    @test flx.ρ == [1, 1, 1]
    @test flx.ρe == [2, 2, 2]
end

# !!! tip
#     A useful pattern for finding tendency definitions
#     is to globally search for, for example,
#     `::Gravity{Momentum` or `::Gravity{`, which
#     (at the time of this writing) yields the results:
#     `function source(s::Gravity{Momentum}, m, args)`

nothing
