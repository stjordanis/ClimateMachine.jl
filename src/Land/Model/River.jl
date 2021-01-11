module River
 
using ..VariableTemplates
using DocStringExtensions
using ..Land
using ..BalanceLaws
import ..BalanceLaws:
    BalanceLaw,
    vars_state,
    flux_first_order!,
    flux_second_order!,
    source!,
    boundary_conditions,
    boundary_state!,
    compute_gradient_argument!,
    compute_gradient_flux!,
    nodal_init_state_auxiliary!,
    init_state_prognostic!,
    nodal_update_auxiliary_state!

export RiverModel, NoRiverModel
#abstract type AbstractRiverModel end

struct NoRiverModel <: BalanceLaw end

struct RiverModel{M} <: BalanceLaw
    mannings::M
end

function RiverModel(mannings)
end

mannings_coeff(river::River Model)

function calculate_velocity(river, x::Real, y::Real, h::Real)
    sx = river.slope_x(x, y)
    sy = river.slope_y(x, y)
    magnitude = h^(2/3) / (river.mannings(x, y)) * sqrt(river.mag_slope(x, y))
    return SVector(sx*magnitude, sy*magnitude, 0)
end

vars_state(water::RiverModel, st::Prognostic, FT) = @vars(height::FT)
vars_state(water::RiverModel, st::Auxiliary, FT) = @vars()
vars_state(water::RiverModel, st::Gradient, FT) = @vars()
vars_state(water::RiverModel, st::GradientFlux, FT) = @vars()

function Land.land_init_aux!(land::LandModel, river::BalanceLaw, aux, geom)
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
    width = river.width(x,y)
    height = state.river.area / width
    v = calculate_velocity(river, x, y ,height)
    Q = state.river.area * v
    flux.river.height = Q
end 


end