module BickleyJet

export BickleyJetModel

using StaticArrays
using ...MPIStateArrays: MPIStateArray
using LinearAlgebra: dot, Diagonal


using ..Ocean
using ...VariableTemplates
using ...Mesh.Geometry
using ...DGMethods
using ...DGMethods.NumericalFluxes
using ...BalanceLaws
using ..Ocean: kinematic_stress, coriolis_parameter

import ...BalanceLaws:
    vars_state,
    init_state_prognostic!,
    init_state_auxiliary!,
    compute_gradient_argument!,
    compute_gradient_flux!,
    flux_first_order!,
    flux_second_order!,
    source!,
    wavespeed,
    boundary_conditions,
    boundary_state!
import ..Ocean:
    ocean_init_state!,
    ocean_init_aux!,
    ocean_boundary_state!,
    _ocean_boundary_state!

using ...Mesh.Geometry: LocalGeometry

×(a::SVector, b::SVector) = StaticArrays.cross(a, b)
⋅(a::SVector, b::SVector) = StaticArrays.dot(a, b)
⊗(a::SVector, b::SVector) = a * b'

abstract type TurbulenceClosure end
struct LinearDrag{T} <: TurbulenceClosure
    λ::L
end
struct ConstantViscosity{T} <: TurbulenceClosure
    ν::T
    κ::T
    function ConstantViscosity{T}(;
        ν = FT(5e3),   # m²/s
        κ = FT(1e3),   # m²/s
    ) where {T <: AbstractFloat}
        return new{T}(ν, κ)
    end
end

abstract type CoriolisForce end
struct fPlaneCoriolis{T} <: CoriolisForce
    fₒ::T
    β::T
    function fPlaneCoriolis{T}(;
        fₒ = T(1e-4), # Hz
        β = T(1e-11), # Hz/m
    ) where {T <: AbstractFloat}
        return new{T}(fₒ, β)
    end
end

abstract type Forcing end
struct KinematicStress{T} <: Forcing
    τₒ::T
    function KinematicStress{T}(; τₒ = T(1e-4)) where {T <: AbstractFloat}
        return new{T}(τₒ)
    end
end

"""
    BickleyJetModel <: BalanceLaw

A `BalanceLaw` for shallow water modeling.

write out the equations here

# Usage

    BickleyJetModel(problem)

"""
struct BickleyJetModel{P, T, A, C, F, FT} <: BalanceLaw
    problem::P
    advection::A
    turbulence::T
    coriolis::C
    forcing::F
    g::FT
    c::FT
    function BickleyJetModel{FT}(
        problem::P,
        turbulence::T,
        advection::A,
        coriolis::C,
        forcing::F,
        g = FT(10), # m/s²
        c = FT(0),
    ) where {FT <: AbstractFloat, P, T, A, C, F}
        return new{P, T, A, C, F, FT}(
            problem,
            turbulence,
            advection,
            coriolis,
            forcing,
            g,
            c,
        )
    end
end
BJModel = BickleyJetModel

function vars_state(m::BJModel, ::Prognostic, T)
    @vars begin
        ρ::T
        ρu::SVector{2, T}
        ρθ::T
    end
end

function init_state_prognostic!(m::BJModel, state::Vars, aux::Vars, localgeo, t)
    ocean_init_state!(m, m.problem, state, aux, localgeo, t)
end

function vars_state(m::BJModel, ::Auxiliary, T)
    @vars begin
        coords::T
    end
end

function init_state_auxiliary!(
    m::BJModel,
    state_auxiliary::MPIStateArray,
    grid,
    direction,
)
    init_state_auxiliary!(
        m,
        (m, A, tmp, geom) -> ocean_init_aux!(m, m.problem, A, geom),
        state_auxiliary,
        grid,
        direction,
    )
end

function vars_state(m::BJModel, ::Gradient, T)
    @vars begin
        ∇u::SVector{2, T}
        ∇θ::T
    end
end

function compute_gradient_argument!(
    model::BJModel,
    grad::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    compute_gradient_argument!(model.turbulence, grad, state, aux, t)
end

compute_gradient_argument!(::LinearDrag, _...) = nothing

@inline function compute_gradient_argument!(
    ::ConstantViscosity,
    grad::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    ρ = state.ρ
    ρu = state.ρu
    ρθ = state.ρθ

    u = ρu / ρ
    θ = ρθ / ρ

    grad.∇u = u
    grad.∇θ = θ

    return nothing
end

function vars_state(m::BJModel, ::GradientFlux, T)
    @vars begin
        ν∇u::SMatrix{3, 2, T, 6}
        κ∇θ::SVector{3, T}
    end
end

function compute_gradient_flux!(
    model::BJModel,
    gradflux::Vars,
    grad::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    compute_gradient_flux!(
        model,
        model.turbulence,
        gradflux,
        grad,
        state,
        aux,
        t,
    )
end

compute_gradient_flux!(::BJModel, ::LinearDrag, _...) = nothing

@inline function compute_gradient_flux!(
    ::BJModel,
    turb::ConstantViscosity,
    gradflux::Vars,
    grad::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    ν = Diagonal(@SVector [turb.ν, turb.ν, -0])
    κ = Diagonal(@SVector [turb.κ, turb.κ, -0])

    ∇u = grad.∇u
    ∇θ = grad.∇θ

    ν∇u = gradflux.ν∇u
    κ∇θ = gradflux.κ∇θ

    ν∇u = -ν * ∇u
    κ∇θ = -κ * ∇θ

    return nothing
end

@inline function flux_first_order!(
    model::BJModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
    direction,
)

    ρ = state.ρ
    ρu = @SVector [state.ρu[1], state.ρu[2], -0]
    ρθ = state.ρθ

    ρₜ = flux.ρ
    ρuₜ = flux.ρu
    θₜ = flux.ρθ

    g = model.g
    h = model.problem.h

    Iʰ = @SMatrix [
        1 -0
        -0 1
        -0 -0
    ]

    ρₜ += ρu
    ρuₜ += 1 // 2 * g * h^2 * Iʰ

    advective_flux!(model, model.advection, flux, state, aux, t)

    return nothing
end

advective_flux!(::BJModel, ::Nothing, _...) = nothing

@inline function advective_flux!(
    ::BJModel,
    ::NonLinearAdvectionTerm,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    ρ = state.ρ
    ρu = @SVector [state.ρu[1], state.ρu[2], -0]
    ρθ = state.ρθ

    ρuₜ = flux.ρu
    ρθₜ = flux.ρθ

    ρuₜ += ρu ⊗ ρu / ρ
    ρθₜ += ρu ⊗ ρθ / ρ

    return nothing
end

function flux_second_order!(
    model::BJModel,
    flux::Grad,
    state::Vars,
    gradflux::Vars,
    ::Vars,
    aux::Vars,
    t::Real,
)
    flux_second_order!(model, model.turbulence, flux, state, gradflux, aux, t)
end

flux_second_order!(::BJModel, ::LinearDrag, _...) = nothing

@inline function flux_second_order!(
    ::BJModel,
    ::ConstantViscosity,
    flux::Grad,
    state::Vars,
    gradflux::Vars,
    aux::Vars,
    t::Real,
)
    ρuₜ = flux.ρu
    ρθₜ = flux.ρθ

    ν∇u = gradflux.ν∇u
    κ∇θ = gradflux.κ∇θ

    ρuₜ += ν∇u
    ρθₜ += κ∇θ

    return nothing
end

@inline function source!(
    model::BJModel,
    source::Vars,
    state::Vars,
    gradflux::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    coriolis_force!(model, model.coriolis, source, state, aux, t)
    forcing_term!(model, model.forcing, source, state, aux, t)
    linear_drag!(model, model.turbulence, source, state, aux, t)

    return nothing
end

coriolis_force!(::BJModel, ::fPlaneCoriolis, _...) = nothing

@inline function coriolis_force!(
    model::BJModel,
    coriolis::fPlaneCoriolis,
    source,
    state,
    aux,
    t,
)
    ρu = @SVector [state.ρu[1], state.ρu[2], -0]
    ρuₜ = source.ρu

    # f × u
    f = [-0, -0, coriolis_parameter(model, coriolis, aux.coords)]
    id = @SVector(1, 2)
    fxρu = (f × ρu)[id]

    ρuₜ -= fxρu

    return nothing
end

forcing_term!(::BJModel, ::Nothing, _...) = nothing

@inline function forcing_term!(
    model::BJModel,
    forcing::KinematicStress,
    source,
    state,
    aux,
    t,
)
    source.ρu += kinematic_stress(model, forcing, aux.coords)

    return nothing
end

linear_drag!(::BJModel, ::ConstantViscosity, _...) = nothing

@inline function linear_drag!(
    ::BJModel,
    turb::LinearDrag,
    source,
    state,
    aux,
    t,
)
    λ = turb.λ
    ρu = state.ρu
    ρuₜ = source.ρu

    ρuₜ -= λ * ρu

    return nothing
end

@inline wavespeed(::BJModel, _...) = m.c

boundary_conditions(model::BJModel) = model.problem.boundary_conditions

"""
    boundary_state!(nf, ::BJModel, args...)

applies boundary conditions for the hyperbolic fluxes
dispatches to a function in OceanBoundaryConditions.jl based on bytype defined by a problem such as SimpleBoxProblem.jl
"""
@inline function boundary_state!(nf, bc, model::BJModel, args...)
    return _ocean_boundary_state!(nf, bc, model, args...)
end

"""
    ocean_boundary_state!(nf, bc::OceanBC, ::BJModel)

splits boundary condition application into velocity
"""
@inline function ocean_boundary_state!(nf, bc::OceanBC, m::BJModel, args...)
    return ocean_boundary_state!(nf, bc.velocity, m, m.turbulence, args...)
    return ocean_boundary_state!(nf, bc.temperature, m, args...)
end

include("bc_velocity.jl")

end
