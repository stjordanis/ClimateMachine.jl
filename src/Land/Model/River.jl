module River
 
using DocStringExtensions
using ..Land
using ..VariableTemplates
import ..BalanceLaws:
    BalanceLaw,
    vars_state,
    flux_first_order!, 
    Prognostic,
    Auxiliary, 
    Gradient,
    GradientFlux
using ...DGMethods: LocalGeometry
using StaticArrays: SVector

export RiverModel, NoRiverModel, river_boundary_flux!, river_boundary_state!

struct NoRiverModel <: BalanceLaw end

struct RiverModel{M,Sx,Sy,MS,W} <: BalanceLaw
    mannings::M
    slope_x::Sx
    slope_y::Sy
    mag_slope::MS
    width::W
end

function RiverModel(
    slope_x::Function,
    slope_y::Function,
    mag_slope::Function,
    width::Function;
    mannings::Function = (x, y) -> convert(eltype(x), 0.03))
    args = (
        slope_x,
        slope_y,
        mag_slope,
        width,
        mannings
    )
    return RiverModel{typeof.(args)...}(args...)
end

function calculate_velocity(river, x::Real, y::Real, h::Real)
    FT = eltype(h)
    sx = FT(river.slope_x(x, y))
    sy = FT(river.slope_y(x, y))
    magnitude = h^FT(2/3) / (river.mannings(x, y) * sqrt(river.mag_slope(x, y)))
    return SVector(sx * magnitude, sy * magnitude, zero(FT))
end

vars_state(river::RiverModel, st::Prognostic, FT) = @vars(area::FT)
vars_state(river::RiverModel, st::Auxiliary, FT) = @vars()
vars_state(river::RiverModel, st::Gradient, FT) = @vars()
vars_state(river::RiverModel, st::GradientFlux, FT) = @vars()

function Land.land_init_aux!(land::LandModel, river::BalanceLaw, aux, geom::LocalGeometry)
end

function Land.compute_gradient_argument!(land::LandModel, river::BalanceLaw, transform::Grad, state, aux, t)
 #   v = calculate_velocity(river, state.river.height,river)
 #   transform.river.grad_Q = state.river.height*v#*river.width(aux.x,aux.y)
end


function Land.land_nodal_update_auxiliary_state!(land::LandModel, river::BalanceLaw, state, aux, t)
end


function flux_first_order!(land::LandModel, river::BalanceLaw, flux::Grad, state::Vars, aux::Vars, t::Real, directions) 
    x = aux.x
    y = aux.y
    width = river.width(x, y)
    height = state.river.area / width
    v = calculate_velocity(river, x, y, height)
    Q = state.river.area * v
    flux.river.area = Q
end 

# boundry conditions 

# General case - to be used with bc::NoBC
function river_boundary_flux!(
    nf,
    bc::Land.AbstractBoundaryConditions,
    m,
    land::LandModel,
    _...,
)
end

function river_boundary_state!(
    nf,
    bc::Land.AbstractBoundaryConditions,
    m,
    land::LandModel,
    _...,
)
end

# Dirichlet BC for River
function river_boundary_flux!(
    nf,
    bc::Land.Dirichlet,
    model::RiverModel,
    land::LandModel,
    state⁺::Vars,
    aux⁺::Vars,
    nM,
    state⁻::Vars,
    aux⁻::Vars,
    t,
    _...,
)
end

function river_boundary_state!(
    nf,
    bc::Land.Dirichlet,
    model::RiverModel,
    land::LandModel,
    state⁺::Vars,
    aux⁺::Vars,
    nM,
    state⁻::Vars,
    aux⁻::Vars,
    t,
    _...,
)
    bc_function = bc.state_bc
    state⁺.river.area = bc_function(aux⁻, t)
end

end