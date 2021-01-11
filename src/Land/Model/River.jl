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
struct RiverModel <: BalanceLaw end

vars_state(water::RiverModel, st::Prognostic, FT) = @vars()
vars_state(water::RiverModel, st::Auxiliary, FT) = @vars()
vars_state(water::RiverModel, st::Gradient, FT) = @vars()
vars_state(water::RiverModel, st::GradientFlux, FT) = @vars()


function Land.land_init_aux!(land, river::BalanceLaw, aux, geom)
end

function Land.compute_gradient_argument!(land, river::BalanceLaw, transform, state, aux, t)
end

function Land.land_nodal_update_auxiliary_state!(land, river::BalanceLaw, state, aux, t)
end




end