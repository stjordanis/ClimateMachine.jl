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
Base.@kwdef struct SoilModelMoisture{AbstractParameterSet, AbstractWater, Fκ, Fiν, Fsν, Fih} <: BalanceLaw #FiS_l, Fiψ_m, Fiψ, Fih
    "Parameters"
    param_set::AbstractParameterSet = param_set
    WF::AbstractWater = waterfunctions()
    #(0.001/(60*60*24)) [m/s] typical value taken from Land Surface Model CLiMA, table 2.2, =0.1cm/day (0.34*1.175e6/(1.175+abs.(aux.h)^4.74))
    initialκ::Fκ       = (aux) -> K_sat
  # Define initial and boundary condition parameters
    initialν::Fiν = (state, aux) -> ν_0 # [m3/m3] constant water content in soil
    surfaceν::Fsν = (state, aux, t) -> ν_surface #267 # [m3/m3] constant flux at surface

  # Define initial and boundary condition parameters
 # initialS_l::FiS_l = (aux) -> 0.1
 # initialψ_m::Fiψ_m = (aux) -> -1
 # initialψ::Fiψ = (aux) ->  -1
 # [m3/m3] constant water content in soil
    initialh::Fih = (aux) -> aux.z + ψ_0 # [m3/m3] constant water content in soil

  #surfaceh::Fsh = (state, aux, t) -> 100  #267 # [m3/m3] constant flux at surface

  # Define initial and boundary condition parameters
  #surfaceψ::Fsψ = (state, aux, t) -> 100  #267 # [m3/m3] constant flux at surface

  # Define initial and boundary condition parameters
  #initialθi::Fiθi = (aux, t) -> 0.0 # [m3/m3] constant water content in soil
  #surfaceψ::Fsψ = (state, aux, t) -> 100  #267 # [m3/m3] constant flux at surface

end

# --------------------------------- 3) Define CliMA vars ---------------------------------------

vars_state_auxiliary(::SoilModelMoisture, FT) = @vars(z::FT, h::FT,κ::FT) #θl::FT, S_l::FT, ψ_m::FT, ψ::FT,
vars_state_conservative(::SoilModelMoisture, FT) = @vars(ν::FT) #, θi::FT)
vars_state_gradient(::SoilModelMoisture, FT) = @vars(h::FT)
vars_state_gradient_flux(::SoilModelMoisture, FT) = @vars(∇h::SVector{3,FT})

# --------------------------------- 4) CliMA functions needed for simulation -------------------
# ---------------- 4a) Initialization

# Initialize z-Profile ### what role does this play? when?
function init_state_auxiliary!(m::SoilModelMoisture, aux::Vars, geom::LocalGeometry)
  aux.z = geom.coord[3]
  #aux.S_l = m.initialS_l(aux) 
  #aux.ψ_m = m.initialψ_m(aux)x
    #aux.ψ = m.initialψ(aux)
    aux.h = m.initialh(aux) #(aux, 0) #aux.z+m.initialθ(state, aux, t) #^(-1/0.378
    aux.κ = m.initialκ(aux)

  # aux.θl = 
end

# Initialize State variables
function init_state_conservative!(m::SoilModelMoisture, state::Vars, aux::Vars, coords, t::Real)
  state.ν = m.initialν(state, aux) #
#  state.θi = m.initialθi(aux, 0)
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

    # How much water    #theta_water = state.θ + state.θi

    # Get augmented liquid
    #theta_l = augmented_liquid(porosity,S_s,aux.ψ,state.θ)

    # Get effective saturation
    S_l = effective_saturation(porosity,state.ν)
#    if isnan(S_l) | isinf(S_l)
#        println(t, state.ν)
#    end
    # This function calculates pressure head ψ of a soil
    ψ = pressure_head(m.WF.matric_pot, S_l,porosity,S_s,state.ν)
#    if isnan(ψ) | isinf(ψ)
#        println(t, S_l,S_s,flag,state.ν)
#    end
    # Get hydraulic head
    aux.h = hydraulic_head(aux.z,ψ)
    aux.κ = hydraulic_conductivity(m.WF.hydraulic_cond,K_sat,S_l,hydraulic_head(aux.z,ψ), aux.z)#, "Havercamp")
#    aux.κ = ksat_function(K_sat,hydraulic_head(aux.z,ψ), aux.z)
    #aux.θl = hydraulic_head(aux.z,aux.ψ)
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

    # How much water
    #theta_water = state.θ + state.θi

    # Get augmented liquid
    #theta_l = augmented_liquid(porosity,S_s,aux.ψ,state.θ)

    ## Get effective saturation
    S_l = effective_saturation(porosity,state.ν)

    ## This function calculates pressure head ψ of a soil
    ψ = pressure_head(m.WF.matric_pot, S_l,porosity,S_s,state.ν)
#    if isnan(ψ) | isinf(ψ)
#        println(t, S_l,S_s,flag,state.ν)
#    end
    # Get hydraulic head
    transform.h = hydraulic_head(aux.z,ψ)

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
   flux.ν -= aux.κ * diffusive.∇h
   #if aux.z == 0
   #end
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


# ---------------- 4e) Boundary Conditions

# Boundary condition function
function boundary_state!(nf, m::SoilModelMoisture, state⁺::Vars, aux⁺::Vars,
                         nM, state⁻::Vars, aux⁻::Vars, bctype, t, _...)
  if bctype == 2
    # surface
      state⁺.ν= m.surfaceν(state⁻, aux⁻, t)
  elseif bctype == 1
    # bottom
    nothing
  end
end

# Boundary condition function
function boundary_state!(nf, m::SoilModelMoisture, state⁺::Vars, diff⁺::Vars,
                         aux⁺::Vars, n̂, state⁻::Vars, diff⁻::Vars, aux⁻::Vars,
                         bctype, t, _...)
  if bctype == 2
    # surface
    state⁺.ν = m.surfaceν(state⁻, aux⁻, t)
  elseif bctype == 1
    # bottom
    #nothing
    #diff⁺.∇h = -diff⁻.∇h
    diff⁺.∇h = -n̂*1
  end
end
