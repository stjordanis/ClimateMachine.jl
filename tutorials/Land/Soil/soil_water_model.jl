# --------------------------------- CLIMA SOIL MODEL -----------------------
# CLIMA_SoilWater.jl: This model simulates soil water dynamics for the CliMA model

"""
Soil Water Model

Computes diffusive flux `F` in:

∂y / ∂t = ∇ ⋅ Flux + Source

```
 ∂(θ)    ∂      ∂h
------ = --(k * --)
  ∂t     ∂z     ∂z
```
where
 - `θ` is the volumetric water content of soil (m³/m³), this is state var.
 - `k` is the hydraulic conductivity (m/s)
 - `h` is the hydraulic head or water potential (m), it is a function of θ
 - `z` is the depth (m)

To write this in the form
```
∂Y
-- + ∇⋅F(Y,t) = 0
∂t
```
we write `Y = θ` and `F(Y, t) =-k ∇h`.

"""

# --------------------------------- 1) Import/Export Needed Functions -----------------------

# Add necessary CliMA functions and sub-routines
using StaticArrays
using ClimateMachine.VariableTemplates
import ClimateMachine.DGmethods: BalanceLaw,
                        vars_state_auxiliary,
                        vars_state_conservative,
                        vars_state_gradient,
                        vars_state_gradient_flux,
                        flux_first_order!,
                        flux_second_order!,
                        source!,
                        compute_gradient_argument!,
                        compute_gradient_flux!,
                        update_auxiliary_state!,
                        nodal_update_auxiliary_state!,
                        init_state_auxiliary!,
                        init_state_conservative!,
                        boundary_state!,
                        wavespeed,
                        LocalGeometry


# --------------------------------- 2) Define Structs ---------------------------------------

"""
Introduce needed variables into SoilModel struct

From Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
"""
Base.@kwdef struct SoilModelMoisture{Fκ, Fiθ, Fsθ, Fih, Fiψ, Fiθi} <: BalanceLaw #, Fsh, Fiψ, Fsψ} <: BalanceLaw
  # Define kappa (hydraulic conductivity)
  #(0.001/(60*60*24)) [m/s] typical value taken from Land Surface Model CLiMA, table 2.2, =0.1cm/day (0.34*1.175e6/(1.175+abs.(aux.h)^4.74))
  K_s::Fκ         = (state, aux, t) -> 1e-3#(1e-3*(0.34/(60*60))*1.175e6/((1.175e6+abs.(aux.h-aux.z)^4.74)))

  # Define initial and boundary condition parameters
  initialθ::Fiθ = (aux, t) -> 0.1 # [m3/m3] constant water content in soil
  surfaceθ::Fsθ = (state, aux, t) -> 0.267 #267 # [m3/m3] constant flux at surface

  # Define initial and boundary condition parameters
  initialh::Fih = (aux, t) -> -1 #- aux.z # [m3/m3] constant water content in soil
  #surfaceh::Fsh = (state, aux, t) -> 100  #267 # [m3/m3] constant flux at surface

  # Define initial and boundary condition parameters
  initialψ::Fiψ = (aux, t) -> -1 - aux.z # [m3/m3] constant water content in soil
  #surfaceψ::Fsψ = (state, aux, t) -> 100  #267 # [m3/m3] constant flux at surface

  # Define initial and boundary condition parameters
  initialθi::Fiθi = (aux, t) -> 0.0 # [m3/m3] constant water content in soil
  #surfaceψ::Fsψ = (state, aux, t) -> 100  #267 # [m3/m3] constant flux at surface

end

# --------------------------------- 3) Define CliMA vars ---------------------------------------

vars_state_auxiliary(::SoilModelMoisture, FT) = @vars(z::FT, h::FT , ψ::FT)
vars_state_conservative(::SoilModelMoisture, FT) = @vars(θ::FT, θi::FT)
vars_state_gradient(::SoilModelMoisture, FT) = @vars(h::FT)
vars_state_gradient_flux(::SoilModelMoisture, FT) = @vars(∇h::SVector{3,FT})

# --------------------------------- 4) CliMA functions needed for simulation -------------------
# ---------------- 4a) Update states

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
  #nodal_update_auxiliary_state
  m::SoilModelMoisture,
  state::Vars,
  aux::Vars,
  t::Real)

    # How much water
    #theta_water = state.θ + state.θi

    # Get augmented liquid
    theta_l = augmented_liquid(porosity,S_s,aux.ψ,theta_water)

    # Get effective saturation
    S_l = effective_saturation(porosity,theta_l)   # 0.2

    # Get matric potential
    ψ_m = matric_potential(flag,S_l)

    # This function calculates pressure head ψ of a soil
    aux.ψ = pressure_head(ψ_m,S_l,porosity,S_s,theta_l)

    # Get hydraulic head
    aux.h = hydraulic_head(aux.z,aux.ψ)
    # transform.h = aux.z+((-1/2.7)*(state.θ/1.987)^(-1/3.96))*(1-state.θ/1.987)^(1/3.96)

end

# ---------------- 4b) Calculate state and derivative of theta

# Calculate h based on state variable
function compute_gradient_argument!(
    m::SoilModelMoisture,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)

    # How much water
    theta_water = state.θ + state.θi

    # Get augmented liquid
    theta_l = augmented_liquid(porosity,S_s,aux.ψ,theta_water)

    # Get effective saturation
    S_l = effective_saturation(porosity,theta_l)   # 0.2

    # Get matric potential
    ψ_m = matric_potential(flag,S_l)

    # This function calculates pressure head ψ of a soil
    aux.ψ = pressure_head(ψ_m,S_l,porosity,S_s,theta_l)

    # Get hydraulic head
    transform.h = hydraulic_head(aux.z,aux.ψ)

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
  diffusive.∇h = ∇transform.h
end

# Calculate thermal flux (non-diffusive)
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
  # Flux of water
   flux.θ -= m.K_s(state, aux, t) * diffusive.∇h
   if aux.z == 0
   end
end

# ---------------- 4c) Extra Sources
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

## Update sources for ice and liquid
#if state.θi > 0
#  # Convert liquid water to ice by freezing (or vice versa)
#  F_T = calculate_frozen_water(state.θ,state.θi,soil_T)
#else
#  F_T = 0;
#end
#
## Source of ice
#source.θi = F_T/917 # rho_i = 0.917 # g cm-3, density of ice
#source.θ = -F_T/997 # rho_l = 0.997 # g cm-3, density of water

end

# ---------------- 4d) Initialization

# Initialize z-Profile ### what role does this play? when?
function init_state_auxiliary!(m::SoilModelMoisture, aux::Vars, geom::LocalGeometry)
  aux.z = geom.coord[3]
  aux.h = m.initialh(aux, 0) #aux.z+m.initialθ(state, aux, t) #^(-1/0.378))*(-0.3020)
  aux.ψ = m.initialψ(aux, 0)
end

# Initialize State variables from T to internal energy
function init_state_conservative!(m::SoilModelMoisture, state::Vars, aux::Vars, coords, t::Real)
  state.θ = m.initialθ(aux, 0)
#  state.θi = m.initialθi(aux, 0)
end

# ---------------- 4e) Boundary Conditions

# Boundary condition function
function boundary_state!(nf, m::SoilModelMoisture, state⁺::Vars, aux⁺::Vars,
                         nM, state⁻::Vars, aux⁻::Vars, bctype, t, _...)
  if bctype == 1
    # surface
    state⁺.θ= m.surfaceθ(state⁻, aux⁻, t)
  elseif bctype == 2
    # bottom
    nothing
  end
end

# Boundary condition function
function boundary_state!(nf, m::SoilModelMoisture, state⁺::Vars, diff⁺::Vars,
                         aux⁺::Vars, nM, state⁻::Vars, diff⁻::Vars, aux⁻::Vars,
                         bctype, t, _...)
  if bctype == 1
    # surface
    state⁺.θ = m.surfaceθ(state⁻, aux⁻, t)
  elseif bctype == 2
    # bottom
    #nothing
    diff⁺.∇h = -diff⁻.∇h
    #diff⁺.∇θ = -diff⁻.∇θ
  end
end
