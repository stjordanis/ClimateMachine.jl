module Runoff

using LinearAlgebra
using DocStringExtensions


using ...VariableTemplates
using ...Land: SoilModel

export AbstractPrecipModel,
    DrivenConstantPrecip,
    AbstractSurfaceRunoffModel,
    NoRunoff,
    compute_surface_state_bc,
    compute_surface_grad_bc,
    CoarseGridRunoff

"""
    AbstractPrecipModel{FT <: AbstractFloat}
"""
abstract type AbstractPrecipModel{FT <: AbstractFloat} end

"""
    DrivenConstantPrecip{FT, F} <: AbstractPrecipModel{FT}

Instance of a precipitation distribution where the precipication value
is constant across the domain. However, this value can change in time.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct DrivenConstantPrecip{FT, F} <: AbstractPrecipModel{FT}
    "Mean precipitation in grid"
    mp::F
    function DrivenConstantPrecip{FT}(mp::F) where {FT, F}
        new{FT, F}(mp)
    end
end

function (dcp::DrivenConstantPrecip{FT})(t::Real) where {FT}
    return FT(dcp.mp(t))
end

"""
    AbstractSurfaceRunoffModel

Abstract type for different surface runoff models. Currently, only
`NoRunoff` is supported.
"""
abstract type AbstractSurfaceRunoffModel end

"""
    NoRunoff <: AbstractSurfaceRunoffModel

Chosen when no runoff is to be modeled.
"""
struct NoRunoff <: AbstractSurfaceRunoffModel end

"""
    CoarseGridRunoff <: AbstractSurfaceRunoffModel

Chosen when no subgrid effects are to be modeled.
"""
struct CoarseGridRunoff <: AbstractSurfaceRunoffModel end


"""
    function compute_surface_grad_bc(soil::SoilModel,
                                     runoff_model::CoarseGridRunoff,
                                     precip_model::AbstractPrecipModel,
                                     n̂,
                                     state⁻::Vars,
                                     diff⁻::Vars,
                                     t::Real
                                     )

Given a runoff model and a precipitation distribution function, compute 
the surface water Neumann BC. This can be a function of time, and state.
"""
function compute_surface_grad_bc(
    soil::SoilModel,
    runoff_model::CoarseGridRunoff,
    precip_model::AbstractPrecipModel,
    n̂,
    state⁻::Vars,
    diff⁻::Vars,
    t::Real,
)
    FT = eltype(state⁻)
    incident_water_flux = precip_model(t)
    if incident_water_flux < -norm(diff⁻.soil.water.K∇h) # More negative if both are negative, 
        #ponding BC
        K∇h⁺ = diff⁻.soil.water.K∇h
    else
        K∇h⁺ = n̂ * (-FT(2) * incident_water_flux) - diff⁻.soil.water.K∇h
    end
    return K∇h⁺
end



"""
    function compute_surface_grad_bc(soil::SoilModel,
                                     runoff_model::NoRunoff,
                                     precip_model::AbstractPrecipModel,
                                     n̂,
                                     state⁻::Vars,
                                     diff⁻::Vars,
                                     t::Real
                                     )

Given a runoff model and a precipitation distribution function, compute 
the surface water Neumann BC. This can be a function of time, and state.
"""
function compute_surface_grad_bc(
    soil::SoilModel,
    runoff_model::NoRunoff,
    precip_model::AbstractPrecipModel,
    n̂,
    state⁻::Vars,
    diff⁻::Vars,
    t::Real,
)
    FT = eltype(state⁻)
    incident_water_flux = precip_model(t)
    K∇h⁺ = n̂ * (-FT(2) * incident_water_flux) - diff⁻.soil.water.K∇h
    return K∇h⁺
end




"""
    function compute_surface_state_bc(soil::SoilModel,
                                      runoff_model::CoarseGridRunoff,
                                      state⁻::Vars,
                                      )

Given a runoff model and a precipitation distribution function, compute 
the surface water Dirichlet BC. This can be a function of time, and state.
"""
function compute_surface_state_bc(
    soil::SoilModel,
    runoff_model::CoarseGridRunoff,
    state⁻::Vars,
)
    FT = eltype(state⁻)
    bc_value = soil.param_functions.porosity - state⁻.soil.water.θ_i
    ϑ_l⁺ = bc_value
    return ϑ_l⁺
end



"""
    function compute_surface_state_bc(soil::SoilModel,
                                      runoff_model::NoRunoff,
                                      state⁻::Vars,
                                      )

Given a runoff model and a precipitation distribution function, compute 
the surface water Dirichlet BC. This can be a function of time, and state.
"""
function compute_surface_state_bc(
    soil::SoilModel,
    runoff_model::NoRunoff,
    state⁻::Vars,
)
    return state⁻.soil.water.ϑ_l
end

end
