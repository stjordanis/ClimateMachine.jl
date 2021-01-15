"""
    ocean_boundary_state!(::NumericalFluxFirstOrder, ::Impenetrable{FreeSlip}, ::BJModel)

apply free slip boundary condition for velocity
sets reflective ghost point
"""
@inline function ocean_boundary_state!(
    ::NumericalFluxFirstOrder,
    ::Impenetrable{FreeSlip},
    ::BJModel,
    ::TurbulenceClosure,
    state⁺,
    aux⁺,
    n⁻,
    state⁻,
    aux⁻,
    t,
    args...,
)
    state⁺.ρ = state⁻.ρ

    ρu⁻ = @SVector [state⁻.ρu[1], state⁻.ρu[2], -0]
    ρu⁺ = ρu⁻ - 2 * n⁻ ⋅ ρu⁻ .* SVector(n⁻)

    id = @SVector(1, 2)
    state⁺.ρu = ρu⁺[id]

    return nothing
end

"""
    ocean_boundary_state!(::Union{NumericalFluxGradient, NumericalFluxSecondOrder}, ::Impenetrable{FreeSlip}, ::BJModel)

no second order flux computed for linear drag
"""
ocean_boundary_state!(
    ::Union{NumericalFluxGradient, NumericalFluxSecondOrder},
    ::VelocityBC,
    ::BJModel,
    ::LinearDrag,
    _...,
) = nothing

"""
    ocean_boundary_state!(::NumericalFluxGradient, ::Impenetrable{FreeSlip}, ::BJModel)

apply free slip boundary condition for velocity
sets non-reflective ghost point
"""
function ocean_boundary_state!(
    ::NumericalFluxGradient,
    ::Impenetrable{FreeSlip},
    ::BJModel,
    ::ConstantViscosity,
    state⁺,
    aux⁺,
    n⁻,
    state⁻,
    aux⁻,
    t,
    args...,
)
    state⁺.ρ = state⁻.ρ

    ρu⁻ = @SVector [state⁻.ρu[1], state⁻.ρu[2], -0]
    ρu⁺ = ρu⁻ - n⁻ ⋅ ρu⁻ .* SVector(n⁻)

    id = @SVector(1, 2)
    state⁺.ρu = ρu⁺[id]

    return nothing
end

"""
    shallow_normal_boundary_flux_second_order!(::NumericalFluxSecondOrder, ::Impenetrable{FreeSlip}, ::BJModel)

apply free slip boundary condition for velocity
apply zero numerical flux in the normal direction
"""
function ocean_boundary_state!(
    ::NumericalFluxSecondOrder,
    ::Impenetrable{FreeSlip},
    ::BJModel,
    ::ConstantViscosity,
    state⁺,
    gradflux⁺,
    aux⁺,
    n⁻,
    state⁻,
    gradflux⁻,
    aux⁻,
    t,
    args...,
)
    state⁺.ρu = state⁻.ρu
    gradflux⁺.ν∇u = n⁻ * (@SVector [-0, -0])'

    return nothing
end

"""
    ocean_boundary_state!(::NumericalFluxFirstOrder, ::Impenetrable{NoSlip}, ::BJModel)

apply no slip boundary condition for velocity
sets reflective ghost point
"""
@inline function ocean_boundary_state!(
    ::NumericalFluxFirstOrder,
    ::Impenetrable{NoSlip},
    ::BJModel,
    ::TurbulenceClosure,
    state⁺,
    aux⁺,
    n⁻,
    state⁻,
    aux⁻,
    t,
    args...,
)
    state⁺.ρ = state⁻.ρ
    state⁺.ρu = -state⁻.ρu

    return nothing
end

"""
    ocean_boundary_state!(::NumericalFluxGradient, ::Impenetrable{NoSlip}, ::BJModel)

apply no slip boundary condition for velocity
set numerical flux to zero for U
"""
@inline function ocean_boundary_state!(
    ::NumericalFluxGradient,
    ::Impenetrable{NoSlip},
    ::BJModel,
    ::ConstantViscosity,
    state⁺,
    aux⁺,
    n⁻,
    state⁻,
    aux⁻,
    t,
    args...,
)
    FT = eltype(state⁺)
    state⁺.ρu = @SVector zeros(FT, 2)

    return nothing
end

"""
    ocean_boundary_state!(::NumericalFluxSecondOrder, ::Impenetrable{NoSlip}, ::BJModel)

apply no slip boundary condition for velocity
sets ghost point to have no numerical flux on the boundary for U
"""
@inline function ocean_boundary_state!(
    ::NumericalFluxSecondOrder,
    ::Impenetrable{NoSlip},
    ::BJModel,
    ::ConstantViscosity,
    state⁺,
    gradflux⁺,
    aux⁺,
    n⁻,
    state⁻,
    gradflux⁻,
    aux⁻,
    t,
    args...,
)
    state⁺.ρu = -state⁻.ρu
    gradflux⁺.ν∇u = gradflux⁻.ν∇u

    return nothing
end

"""
    ocean_boundary_state!(::Union{NumericalFluxFirstOrder, NumericalFluxGradient}, ::Penetrable{FreeSlip}, ::BJModel)

no mass boundary condition for penetrable
"""
ocean_boundary_state!(
    ::Union{NumericalFluxFirstOrder, NumericalFluxGradient},
    ::Penetrable{FreeSlip},
    ::BJModel,
    ::ConstantViscosity,
    _...,
) = nothing

"""
    ocean_boundary_state!(::NumericalFluxSecondOrder, ::Penetrable{FreeSlip}, ::BJModel)

apply free slip boundary condition for velocity
apply zero numerical flux in the normal direction
"""
function ocean_boundary_state!(
    ::NumericalFluxSecondOrder,
    ::Penetrable{FreeSlip},
    ::BJModel,
    ::ConstantViscosity,
    state⁺,
    gradflux⁺,
    aux⁺,
    n⁻,
    state⁻,
    gradflux⁻,
    aux⁻,
    t,
    args...,
)
    state⁺.ρu = state⁻.ρu
    gradflux⁺.ν∇u = n⁻ * (@SVector [-0, -0])'

    return nothing
end

"""
    ocean_boundary_state!(::Union{NumericalFluxFirstOrder, NumericalFluxGradient}, ::Impenetrable{KinematicStress}, ::HBModel)

apply kinematic stress boundary condition for velocity
applies free slip conditions for first-order and gradient fluxes
"""
function ocean_boundary_state!(
    nf::Union{NumericalFluxFirstOrder, NumericalFluxGradient},
    ::Impenetrable{<:KinematicStress},
    model::BJModel,
    turb::TurbulenceClosure,
    args...,
)
    return ocean_boundary_state!(
        nf,
        Impenetrable(FreeSlip()),
        model,
        turb,
        args...,
    )
end

"""
    ocean_boundary_state!(::NumericalFluxSecondOrder, ::Impenetrable{KinematicStress}, ::HBModel)

apply kinematic stress boundary condition for velocity
sets ghost point to have specified flux on the boundary for ν∇u
"""
@inline function ocean_boundary_state!(
    ::NumericalFluxSecondOrder,
    ::Impenetrable{<:KinematicStress},
    model::BJModel,
    state⁺,
    gradflux⁺,
    aux⁺,
    n⁻,
    state⁻,
    gradflux⁻,
    aux⁻,
    t,
)
    state⁺.ρu = state⁻.ρu
    gradflux⁺.ν∇u = n⁻ * kinematic_stress(model, aux⁻.coords, state⁻.ρ)'
    # applies windstress for now, will be fixed in a later PR

    return nothing
end

"""
    ocean_boundary_state!(::Union{NumericalFluxFirstOrder, NumericalFluxGradient}, ::Penetrable{KinematicStress}, ::HBModel)

apply kinematic stress boundary condition for velocity
applies free slip conditions for first-order and gradient fluxes
"""
function ocean_boundary_state!(
    nf::Union{NumericalFluxFirstOrder, NumericalFluxGradient},
    ::Penetrable{<:KinematicStress},
    model::BJModel,
    turb::TurbulenceClosure,
    args...,
)
    return ocean_boundary_state!(
        nf,
        Penetrable(FreeSlip()),
        model,
        turb,
        args...,
    )
end

"""
    ocean_boundary_state!(::NumericalFluxSecondOrder, ::Penetrable{KinematicStress}, ::HBModel)

apply kinematic stress boundary condition for velocity
sets ghost point to have specified flux on the boundary for ν∇u
"""
@inline function ocean_boundary_state!(
    ::NumericalFluxSecondOrder,
    ::Penetrable{<:KinematicStress},
    shallow::BJModel,
    ::TurbulenceClosure,
    state⁺,
    gradflux⁺,
    aux⁺,
    n⁻,
    state⁻,
    gradflux⁻,
    aux⁻,
    t,
)
    state⁺.ρu = state⁻.ρu
    gradflux⁺.ν∇u = n⁻ * kinematic_stress(model, aux⁻.coords, state⁻.ρ)'
    # applies windstress for now, will be fixed in a later PR

    return nothing
end
