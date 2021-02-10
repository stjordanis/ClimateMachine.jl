if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end

using ClimateMachine.BalanceLaws
using ClimateMachine.VariableTemplates
using StaticArrays, Test

import ClimateMachine.BalanceLaws: prognostic_vars, eq_tends, flux

struct MyBalanceLaw <: BalanceLaw end

struct Mass <: PrognosticVariable end
struct Energy <: PrognosticVariable end

prognostic_vars(::MyBalanceLaw) = (Mass(), Energy());

struct Advection{PV} <: TendencyDef{Flux{FirstOrder}, PV} end
struct Source1{PV} <: TendencyDef{Source, PV} end
struct Source2{PV} <: TendencyDef{Source, PV} end

struct Diffusion{PV <: Energy} <: TendencyDef{Flux{SecondOrder}, PV} end

#! format: off
eq_tends(pv::PV, ::MyBalanceLaw, ::Flux{FirstOrder}) where {PV <: Mass} = (Advection{PV}(),);
eq_tends(pv::PV, ::MyBalanceLaw, ::Flux{FirstOrder}) where {PV <: Energy} = (Advection{PV}(),);
eq_tends(pv::PV, ::MyBalanceLaw, ::Flux{SecondOrder}) where {PV <: Mass} = ();
eq_tends(pv::PV, ::MyBalanceLaw, ::Flux{SecondOrder}) where {PV <: Energy} = (Diffusion{PV}(),);
eq_tends(pv::PV, ::MyBalanceLaw, ::Source) where {PV <: Mass} = (Source1{PV}(), Source2{PV}());
eq_tends(pv::PV, ::MyBalanceLaw, ::Source) where {PV <: Energy} = (Source1{PV}(), Source2{PV}());
#! format: on

bl = MyBalanceLaw()
show_tendencies(bl; table_complete = true)

flux(::Advection{Mass}, bl::MyBalanceLaw, args) = args.state.ρ;
flux(::Advection{Energy}, bl::MyBalanceLaw, args) = args.state.ρe;

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

    # `Σfluxes(eq_tends(Mass(), bl, tend), bl, args)` calls
    # `flux(::Advection{Mass}, ...)` defined above:
    flx.ρ = Σfluxes(eq_tends(Mass(), bl, tend), bl, args) .* vec_pad

    # `Σfluxes(eq_tends(Energy(), bl, tend), bl, args)` calls
    # `flux(::Advection{Energy}, ...)` defined above:
    flx.ρe = Σfluxes(eq_tends(Energy(), bl, tend), bl, args) .* vec_pad
    return nothing
end;

FT = Float64; # float type
aux = (); # auxiliary fields
t = 0.0; # time
direction = nothing; # Direction

state = Vars{@vars(ρ::FT, ρe::FT)}([1, 2]);
flx = Grad{@vars(ρ::FT, ρe::FT)}(zeros(MArray{Tuple{3, 2}, FT}));

flux_first_order!(bl, flx, state, aux, t, direction);

@testset "Test results" begin
    @test flx.ρ == [1, 1, 1]
    @test flx.ρe == [2, 2, 2]
end

nothing

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

