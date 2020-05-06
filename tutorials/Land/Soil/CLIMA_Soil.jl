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
 - `p` is matric potential (m), p=-(θ/v)^(-M) with Brooks and Corey Formulation
 - `v` is porosity, typical value :
 - `M` is paramter, typical value for sand: 1/0.378, from  BRAUN 2004, p. 1118

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
using CLIMA.VariableTemplates
import CLIMA.DGmethods: BalanceLaw,
                        vars_state_auxiliary, vars_state_conservative, vars_state_gradient, vars_state_gradient_flux,
                        flux_first_order!, flux_second_order!, source!,
                        compute_gradient_argument!, compute_gradient_flux!, update_auxiliary_state!, nodal_update_auxiliary_state!,
                        init_state_auxiliary!, init_state_conservative!,
                        boundary_state!, wavespeed, LocalGeometry


# --------------------------------- 2) Define Structs ---------------------------------------


# Introduce needed variables into SoilModel struct
Base.@kwdef struct SoilModels{Fk, Fiθ, Fsθ, Fih, Fiψ, Fiθi, Fρc, Fκ, FiT, Fst} <: BalanceLaw #, Fsh, Fiψ, Fsψ} <: BalanceLaw
  # Define kappa (hydraulic conductivity)
  K_s::Fk         = (state, aux, t) -> 1e-3#(1e-3*(0.34/(60*60))*1.175e6/((1.175e6+abs.(aux.h-aux.z)^4.74))) #(0.001/(60*60*24)) [m/s] typical value taken from Land Surface Model CLiMA, table 2.2, =0.1cm/day (0.34*1.175e6/(1.175+abs.(aux.h)^4.74)) 
  
  # Define initial and boundary condition parameters
  initialθ::Fiθ = (aux, t) -> 0.1 # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
  surfaceθ::Fsθ = (state, aux, t) -> 0.267 #267 # [m3/m3] constant flux at surface, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287

  # Define initial and boundary condition parameters
  initialh::Fih = (aux, t) -> -1 #- aux.z # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
  #surfaceh::Fsh = (state, aux, t) -> 100  #267 # [m3/m3] constant flux at surface, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287

  # Define initial and boundary condition parameters
  initialψ::Fiψ = (aux, t) -> -1 - aux.z # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
  #surfaceψ::Fsψ = (state, aux, t) -> 100  #267 # [m3/m3] constant flux at surface, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287

  # Define initial and boundary condition parameters
  initialθi::Fiθi = (aux, t) -> 0.0 # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
  #surfaceψ::Fsψ = (state, aux, t) -> 100  #267 # [m3/m3] constant flux at surface, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287

  # Define heat capacity. This is an input to the model now.
  ρc::Fρc       = (state, aux, t) -> 2.49e6   # [ Sand: ρc = 2.49e6 J m-3 K-1 ; Clay: ρc = 2.61e6 J m-3 K-1 ]
  
  κ::Fκ         = (state, aux, t) -> 2.42     # [ Sand: λ = 2.42 W m-1 K-1 ; Clay: λ = 1.17 W m-1 K-1 ]
  # Define kappa (thermal conductivity). This is an input to the model now.

  # # Define initial and boundary condition parameters
  initialT::FiT = (aux, t) -> 290 # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
  surfaceT::Fst = (state, aux, t) -> 280  # Surface boundary condition. This is an input to the model now.

end


# --------------------------------- 3) Define CliMA vars ---------------------------------------


# Stored in the aux state are:
#   `coord` coordinate points (needed for BCs)
#   `u` advection velocity
#   `D` Diffusion tensor
# vars_state_auxiliary(::SoilModelMoisture, FT) = @vars(z::FT, h::FT , ψ::FT) # p::Fθ stored in dg.auxstate
# vars_state_conservative(::SoilModelMoisture, FT) = @vars(θ::FT, θi::FT) # stored in Q
# vars_state_gradient(::SoilModelMoisture, FT) = @vars(h::FT) # not stored
# vars_state_gradient_flux(::SoilModelMoisture, FT) = @vars(∇h::SVector{3,FT}) # stored in dg.diffstate

vars_state_auxiliary(::SoilModels, FT) = @vars(z::FT, h::FT , ψ::FT, T::FT) # p::Fθ stored in dg.auxstate
vars_state_conservative(::SoilModels, FT) = @vars(ρcT::FT, θ::FT, θi::FT) # stored in Q
vars_state_gradient(::SoilModels, FT) = @vars(h::FT, T::FT) # not stored
vars_state_gradient_flux(::SoilModels, FT) = @vars(∇h::SVector{3,FT}, ∇T::SVector{3,FT}) # stored in dg.diffstate


# --------------------------------- 4) CliMA functions needed for simulation -------------------

# ---------------- 4a) Update states

# Update all auxiliary variables
function update_auxiliary_state!(
    dg::DGModel,
    m::SoilModels,
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
  m::SoilModels,
  state::Vars,
  aux::Vars,
  t::Real)
    # flag = "van Genuchten" # - "Brooks and Corey"
    
    # Soil Matric potential - "van Genuchten"
    if flag == "van Genuchten"
        alpha = 0.02 # m-1
        n = 5
        M = 1 - 1/n 
    elseif flag == "Brooks and Corey"
    # Soil Matric potential - "Brooks and Corey"
        alpha = 0.02 # m-1
        n = 5
        M = 1 - 1/n 
    end
    
    # How much water
    theta_water = state.θ + state.θi
# @show theta_water

    # Get augmented liquid
    theta_l = augmented_liquid(porosity,S_s,aux.ψ,theta_water) 

    # Get effective saturation
    S_l = effective_saturation(porosity,theta_l)   # 0.2

    # Get matric potential
    ψ_m = matric_potential(flag,alpha,S_l,n,M)

    # This function calculates pressure head ψ of a soil
    aux.ψ = pressure_head(ψ_m,S_l,porosity,S_s,theta_l)  
  
    # Get hydraulic head    
    aux.h = hydraulic_head(aux.z,aux.ψ)        
    # transform.h = aux.z+((-1/2.7)*(state.θ/1.987)^(-1/3.96))*(1-state.θ/1.987)^(1/3.96)

    # Update Temperature
    # @show m.ρc(state, aux, t)
    # @show state.ρcT
    # @show state.θi

   # ρc = m.ρc(state, aux, t)
   # aux.T = temperature_calculator(ρc,state.ρcT,state.θi)
    aux.T = temperature_calculator(m.ρc(state, aux, t),state.ρcT,state.θi)
#@show aux.T
end


# ---------------- 4b) Calculate state and derivative of theta

# Calculate h based on state variable
function compute_gradient_argument!(
    m::SoilModels,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)     
    # flag = "van Genuchten" # - "Brooks and Corey"
    
    # Soil Matric potential - "van Genuchten"
    if flag == "van Genuchten"
        alpha = 0.02 # m-1
        n = 5
        M = 1 - 1/n 
    elseif flag == "Brooks and Corey"
    # Soil Matric potential - "Brooks and Corey"
        alpha = 0.02 # m-1
        n = 5
        M = 1 - 1/n 
    end

    # How much water
    theta_water = state.θ + state.θi
# @show theta_water

    # Get augmented liquid
    theta_l = augmented_liquid(porosity,S_s,aux.ψ,theta_water) 

    # Get effective saturation
    S_l = effective_saturation(porosity,theta_l)   # 0.2
# @show S_l

    # Get matric potential
    ψ_m = matric_potential(flag,alpha,S_l,n,M)

    # This function calculates pressure head ψ of a soil
    aux.ψ = pressure_head(ψ_m,S_l,porosity,S_s,theta_l)  
  
    # Get hydraulic head    
    transform.h = hydraulic_head(aux.z,aux.ψ)        

    # Update Temperature
   # ρc = m.ρc(state, aux, t)
   # transform.T = temperature_calculator(ρc,state.ρcT,state.θi)
   transform.T = temperature_calculator(m.ρc(state, aux, t),state.ρcT,state.θi)
# @show m.ρc(state, aux, t)
# @show state.ρcT
# @show transform.T
end

# Gradient of h calculation
function compute_gradient_flux!(
    m::SoilModels,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
  diffusive.∇h = ∇transform.h
  diffusive.∇T = ∇transform.T
end

# Calculate thermal flux (non-diffusive)
function  flux_first_order!(
    m::SoilModels,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )

end

# Calculate water flux (diffusive)
function flux_second_order!(
    m::SoilModels,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
  )
    # Flux of heat
  flux.ρcT -= m.κ(state, aux, t) * diffusive.∇T
      # @show flux.ρcT
  # Flux of water
  flux.θ -= m.K_s(state, aux, t) * diffusive.∇h
     # @show flux.θ
    # @show flux.θ
   if aux.z == 0
   end
end

# ---------------- 4c) Extra Sources

# Introduce sources of energy (e.g. Metabolic heat from microbes)
function source!(
    m::SoilModels,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)

# @show state.θi
# @show state.θ
# @show aux.T

# Update sources for ice and liquid
# Convert liquid water to ice by freezing (or vice versa)
F_T = calculate_frozen_water(state.θ,state.θi,aux.T)
# @show F_T

# Source of ice
source.θi = F_T/917 # rho_i = 0.917 # g cm-3, density of ice
source.θ = -F_T/997 # rho_l = 0.997 # g cm-3, density of water

# @show source.θi
# @show source.θ

end


# ---------------- 4d) Initialization 


# Initialize z-Profile ### what role does this play? when?
function init_state_auxiliary!(m::SoilModels, aux::Vars, geom::LocalGeometry)
  aux.z = geom.coord[3]
  aux.h = m.initialh(aux, 0) #aux.z+m.initialθ(state, aux, t) #^(-1/0.378))*(-0.3020)
  aux.ψ = m.initialψ(aux, 0)
  aux.T = m.initialT(aux, 0)
end

# Initialize State variables from T to internal energy
function init_state_conservative!(m::SoilModels, state::Vars, aux::Vars, coords, t::Real)
  state.θ = m.initialθ(aux, 0)
  state.θi = m.initialθi(aux, 0)
  state.ρcT = m.ρc(state, aux, 0) * m.initialT(aux, 0) 
  # @show m.ρc(state, aux, t)
  # @show 1
end


# ---------------- 4e) Boundary Conditions


# Boundary condition function
function boundary_state!(nf, m::SoilModels, state⁺::Vars, aux⁺::Vars,
                         nM, state⁻::Vars, aux⁻::Vars, bctype, t, _...)
  if bctype == 1
    # surface
    state⁺.ρcT = m.ρc(state⁻, aux⁻, t) * m.surfaceT(state⁻, aux⁻, t)
    state⁺.θ= m.surfaceθ(state⁻, aux⁻, t)
  elseif bctype == 2
    # bottom
    nothing
  end
end

# Boundary condition function
function boundary_state!(nf, m::SoilModels, state⁺::Vars, diff⁺::Vars,
                         aux⁺::Vars, nM, state⁻::Vars, diff⁻::Vars, aux⁻::Vars,
                         bctype, t, _...)
  if bctype == 1
    # surface
    state⁺.ρcT = m.ρc(state⁻, aux⁻, t) * m.surfaceT(state⁻, aux⁻, t)
    state⁺.θ = m.surfaceθ(state⁻, aux⁻, t)
  elseif bctype == 2
    # bottom
    nothing
    #diff⁺.∇h = -diff⁻.∇h
    #diff⁺.∇θ = -diff⁻.∇θ        
  end
end
