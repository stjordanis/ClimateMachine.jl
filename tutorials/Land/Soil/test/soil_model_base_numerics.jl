# --------------------------------- CLIMA SOIL MODEL -----------------------
# CLIMA_SoilWater.jl: This model simulates soil water dynamics for the CliMA model

"""
Soil Water Model

Computes diffusive flux `F` in:

∂y / ∂t = ∇ ⋅ Flux + Source

```
 ∂(ν)      ∂      ∂h
------ = - --(-k *--)
  ∂t       ∂z     ∂z
```
where
 - `ν` is the augmented volumetric water content of soil (m³/m³), this is state var.
 - `k` is the hydraulic conductivity (m/s)
 - `h` is the hydraulic head or water potential (m), it is a function of ν
 - `z` is the depth (m)

To write this in the form
```
∂Y
-- + ∇⋅F(Y,t) = 0
∂t
```
we write `Y = ν` and `F(Y, t) =-k ∇h`.

"""

# --------------------------------- 1) Import/Export Needed Functions -----------------------
import ClimateMachine.DGMethods:
    vars_state_auxiliary,
    vars_state_conservative,
    vars_state_gradient,
    vars_state_gradient_flux,
    source!,
    flux_second_order!,
    flux_first_order!,
    compute_gradient_argument!,
    compute_gradient_flux!,
    update_auxiliary_state!,
    nodal_update_auxiliary_state!,
    init_state_auxiliary!,
    init_state_conservative!,
    boundary_state!

# --------------------------------- 2) Define Structs ---------------------------------------

"""
Introduce needed variables into SoilModel struct

From Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287

"""
Base.@kwdef struct SoilModelMoisture{AbstractParameterSet, AbstractWater, Fκ, Fiν, Fih, Dbc, Nbc} <: BalanceLaw 
    "Parameters"
    param_set::AbstractParameterSet = param_set
    WF::AbstractWater = waterfunctions()
  # Define initial and boundary condition parameters
    initialκ::Fκ       = (aux) -> K_sat
    initialν::Fiν = (state, aux) -> ν_0
    initialh::Fih = (aux) -> aux.z + ψ_0
    dirichlet_bc::Dbc = nothing
    neumann_bc::Nbc = nothing
    



end

# --------------------------------- 3) Define CliMA vars ---------------------------------------

vars_state_auxiliary(::SoilModelMoisture, FT) = @vars(z::FT, h::FT,κ::FT) #θl::FT, S_l::FT, ψ_m::FT, ψ::FT,
vars_state_conservative(::SoilModelMoisture, FT) = @vars(ν::FT) #, θi::FT)
vars_state_gradient(::SoilModelMoisture, FT) = @vars(h::FT)
vars_state_gradient_flux(::SoilModelMoisture, FT) = @vars(κ∇h::SVector{3,FT})#really, the flux is - κ∇h

# --------------------------------- 4) CliMA functions needed for simulation -------------------
# ---------------- 4a) Initialization

# Initialize z-Profile ### what role does this play? when?
function init_state_auxiliary!(m::SoilModelMoisture, aux::Vars, geom::LocalGeometry)
    aux.z = geom.coord[3]
    aux.h = m.initialh(aux) 
    aux.κ = m.initialκ(aux)

  # aux.θl = 
end

# Initialize State variables
function init_state_conservative!(m::SoilModelMoisture, state::Vars, aux::Vars, coords, t::Real)
  state.ν = m.initialν(state, aux) #
end


# ---------------- 4b) Update states

# Update all auxiliary variables
function update_auxiliary_state!(
    dg::DGModel,
    m::SoilModelMoisture,
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
)
  nodal_update_auxiliary_state!(soil_nodal_update_aux!, dg, m, Q, t, elems)
  return true
end

# Update all auxiliary nodes
function  soil_nodal_update_aux!(
  m::SoilModelMoisture,
  state::Vars,
  aux::Vars,
  t::Real)

    # Get effective saturation
#    # For debugging only:
#    if state.ν < 0
#        println("Negative moisture content at ", aux.z)
#    end
    S_l = effective_saturation(porosity,state.ν)
    # This function calculates pressure head ψ of a soil
    ψ = pressure_head(m.WF.matric_pot,porosity,S_s,state.ν)
    # Get hydraulic head
    aux.h = hydraulic_head(aux.z,ψ)
    aux.κ = hydraulic_conductivity(m.WF.hydraulic_cond,K_sat,S_l,ψ)
end

# ---------------- 4c) Calculate state and derivative of theta

# Calculate h based on state variable
function compute_gradient_argument!(
    m::SoilModelMoisture,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    # Get effective saturation
    S_l = effective_saturation(porosity,state.ν)

    # This function calculates pressure head ψ of a soil
    ψ = pressure_head(m.WF.matric_pot, porosity,S_s,state.ν)
    # Get hydraulic head
    transform.h = hydraulic_head(aux.z,ψ)#This can't be aux.h here. Why?

end

# Gradient of h calculation
function compute_gradient_flux!(
    m::SoilModelMoisture,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
  diffusive.κ∇h = aux.κ*∇transform.h
end

# Calculate water flux (non-diffusive)
function  flux_first_order!(
    m::SoilModelMoisture,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
end

# Calculate water flux (diffusive)
function flux_second_order!(
    m::SoilModelMoisture,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
  )
   flux.ν -= diffusive.κ∇h
end

# ---------------- 4d) Extra Sources
# Introduce sources of energy (e.g. Metabolic heat from microbes)
function source!(
    m::SoilModelMoisture,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)

end


# ---------------- 4e) Boundary Conditions

# Boundary condition function
function boundary_state!(nf, m::SoilModelMoisture, state⁺::Vars, aux⁺::Vars,
                         nM, state⁻::Vars, aux⁻::Vars, bctype, t, _...)

    bc = m.dirichlet_bc
    if bctype == 2
        top_boundary_conditions!(bc, state⁺, aux⁺, state⁻, aux⁻,  t)
    elseif bctype == 1
        bottom_boundary_conditions!(bc, state⁺, aux⁺, state⁻, aux⁻, t)
    end
end

# Boundary condition function
function boundary_state!(nf, m::SoilModelMoisture, state⁺::Vars, diff⁺::Vars,
                         aux⁺::Vars, n̂, state⁻::Vars, diff⁻::Vars, aux⁻::Vars,
                         bctype, t, _...)
    #do I need to redefine the dirichlet BC here? We did that in the original soil model
    bc = m.neumann_bc
    if bctype == 2
        top_boundary_conditions!(bc, state⁺, diff⁺, aux⁺, n̂, state⁻, diff⁻, aux⁻, t)
    elseif bctype == 1
        bottom_boundary_conditions!(bc, state⁺, diff⁺, aux⁺, n̂, state⁻, diff⁻, aux⁻, t)
  end
end


abstract type bc_values{FT <: AbstractFloat} end

mutable struct Dirichlet{FT} <: bc_values{FT}
    surface_state::FT
    bottom_state::FT
end
#scalars only
mutable struct Neumann{FT} <: bc_values{FT}
    "this will be multiplied by the normal vector - TBD is this plus or minus"
    surface_flux::FT
    "this will be multiplied by the normal vector  - I think it points out of the domain; in minus zhat."
    bottom_flux::FT
end

function top_boundary_conditions!(bc::Neumann{FT}, state⁺::Vars, diff⁺::Vars,
                                  aux⁺::Vars, n̂, state⁻::Vars, diff⁻::Vars, aux⁻::Vars,
                                  t) where {FT}
    if !isnan(bc.surface_flux)
        diff⁺.κ∇h = n̂*bc.surface_flux*aux⁻.κ
    else
        nothing
    end
end

function top_boundary_conditions!(bc::Dirichlet{FT}, state⁺::Vars, aux⁺::Vars,
                                  state⁻::Vars, aux⁻::Vars, t
                                  ) where {FT}
    if !isnan(bc.surface_state)
        state⁺.ν = bc.surface_state
    else
        nothing
    end
end


function bottom_boundary_conditions!(bc::Neumann{FT}, state⁺::Vars, diff⁺::Vars,
                                     aux⁺::Vars, n̂, state⁻::Vars, diff⁻::Vars, aux⁻::Vars,t
                                     ) where {FT}
    if !isnan(bc.bottom_flux)
        diff⁺.κ∇h = n̂*bc.bottom_flux*aux⁻.κ
    else
        nothing
    end
end

function bottom_boundary_conditions!(bc::Dirichlet{FT}, state⁺::Vars, aux⁺::Vars,
                                     state⁻::Vars, aux⁻::Vars,t
                                     ) where {FT}
    if !isnan(bc.bottom_state)
        state⁺.ν = bc.bottom_state
    else
        nothing
    end
end
