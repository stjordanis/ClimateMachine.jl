# --------------------------------- CLIMA SOIL MODEL -----------------------
# CLIMA_Soil.jl: This model simulates soil heat + water dynamics for the CliMA model

#= ---------------------------------------------------------------------------------
Soil Heat Model
Computes diffusive flux `F` in:
∂y / ∂t = ∇ ⋅ Flux + Source
```
∂(ρcT)   ∂      ∂T
------ = --(λ * --)
  ∂t     ∂z     ∂z
```
where
 - `ρ` is the density of the soil (kg/m³)
 - `c` is the soil heat capacity (J/(kg K))
 - `λ` is the thermal conductivity (W/(m K))
To write this in the form
```
∂Y
-- + ∇⋅F(Y,t) = 0
∂t
```
we write `Y = ρcT` and `F(Y, t) = -λ ∇T`.
"""
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

=# # ---------------------------------------------------------------------------------

# --------------------------------- 1) Import/Export Needed Functions ----------------

# Add necessary CliMA functions and sub-routines
using StaticArrays
using CLIMA.VariableTemplates
import CLIMA.DGmethods: BalanceLaw,
                        vars_aux, vars_state, vars_gradient, vars_diffusive,
                        flux_nondiffusive!, flux_diffusive!, source!,
                        gradvariables!, diffusive!, update_aux!, nodal_update_aux!,
                        init_aux!, init_state!,
                        boundary_state!, wavespeed, LocalGeometry


# --------------------------------- 2) Define Structs ---------------------------------------

# Introduce needed variables into SoilModelFull struct
Base.@kwdef struct SoilModelFull{Fρc, Fκ, FiT, Fst, Fk, Fiθ, Fsθ, Fih, Fiψ, Fiθi} <: BalanceLaw
 # HEAT
  # Define heat capacity. This is an input to the model now.
  ρc::Fρc       = (state, aux, t) -> 2.49e6   # [ Sand: ρc = 2.49e6 J m-3 K-1 ; Clay: ρc = 2.61e6 J m-3 K-1 ]
  # Replace this with a function that calculates heat capacity (based on liquid+ice)
  # OR Replace this with tabulated values of heat capacity (based on liquid+ice)

  # Define kappa (thermal conductivity). This is an input to the model now.
  κ::Fκ         = (state, aux, t) -> 2.42     # [ Sand: λ = 2.42 W m-1 K-1 ; Clay: λ = 1.17 W m-1 K-1 ]

  # Define initial and boundary condition parameters
  initialT::FiT = (aux, t) -> 273.15 + 2.0 # Initial Temperature. This is an input to the model now.
  surfaceT::Fst = (state, aux, t) -> (273.15 + 2.0) # Surface boundary condition. This is an input to the model now.

# WATER
  # Define kappa (hydraulic conductivity)
  K_s::Fκ       = (state, aux, t) -> 1e-3#(1e-3*(0.34/(60*60))*1.175e6/((1.175e6+abs.(aux.h-aux.z)^4.74))) #(0.001/(60*60*24)) [m/s] typical value taken from Land Surface Model CLiMA, table 2.2, =0.1cm/day (0.34*1.175e6/(1.175+abs.(aux.h)^4.74)) 
  
  # Define initial and boundary condition parameters
  initialθ::Fiθ = (aux, t) -> 0.1 # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
  surfaceθ::Fsθ = (state, aux, t) -> 0.267 #267 # [m3/m3] constant flux at surface, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287

  # Define initial and boundary condition parameters
  initialh::Fih = (aux, t) -> -1 #- aux.z # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
 
  # Define initial and boundary condition parameters
  initialψ::Fiψ = (aux, t) -> -1 - aux.z # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
 
  # Define initial and boundary condition parameters
  initialθi::Fiθi = (aux, t) -> 0.0 # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
 
end


# --------------------------------- 3) Define CliMA vars ---------------------------------------

# Stored in the aux state are:
#   `coord` coordinate points (needed for BCs)
#   `u` advection velocity
#   `D` Diffusion tensor

vars_aux(::SoilModelFull, FT) = @vars(z::FT, T::FT, h::FT , ψ::FT) # stored dg.auxstate
vars_state(::SoilModelFull, FT) = @vars(ρcT::FT, θ::FT, θi::FT) # stored in Q
vars_gradient(::SoilModelFull, FT) = @vars(T::FT, ψ::FT) # not stored
vars_diffusive(::SoilModelFull, FT) = @vars(∇T::SVector{3,FT}, ∇h::SVector{3,FT}) # stored in dg.diffstate

# vars_aux(::SoilModelFull, FT) = @vars(z::FT, T::FT) # stored dg.auxstate
# vars_state(::SoilModelFull, FT) = @vars(ρcT::FT, θ::FT, θi::FT) # stored in Q
# vars_gradient(::SoilModelFull, FT) = @vars(T::FT) # not stored
# vars_diffusive(::SoilModelFull, FT) = @vars(∇T::SVector{3,FT}) # stored in dg.diffstate

# vars_aux(::SoilModelFull, Fθ) = @vars(z::Fθ, h::Fθ , ψ::Fθ) # p::Fθ stored in dg.auxstate
# vars_state(::SoilModelFull, Fθ) = @vars(θ::Fθ, θi::Fθ) # stored in Q
# vars_gradient(::SoilModelFull, Fθ) = @vars(h::Fθ) # not stored
# vars_diffusive(::SoilModelFull, Fθ) = @vars(∇h::SVector{3,Fθ}) # stored in dg.diffstate


# --------------------------------- 4) CliMA functions needed for simulation -------------------

# ---------------- 4a) Update states

# Update all auxiliary variables
function update_aux!(
    dg::DGModel,
    m::SoilModelFull,
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
)
  nodal_update_aux!(soil_nodal_update_aux!, dg, m, Q, t, elems)
  return true
end
# Update all auxiliary nodes
function soil_nodal_update_aux!(
  m::SoilModelFull,
  state::Vars,
  aux::Vars,
  t::Real)

# HEAT
  #aux.T = state.ρcT / m.ρc(state, aux, t) # TODO: figure out why can't use aux.T here
    # aux.T = temperature_calculator(m.ρc(state, aux, t),state.ρcT,state.θi)
  aux.T = temperature_calculator(m.ρc(state, aux, t),state.ρcT,state.θi)

# WATER
   # Soil Matric potential - "van Genuchten"
    if flag == "van Genuchten"
        alpha = 2 # m-1
        n = 5
        m = 1 - 1/n 
    elseif flag == "Brooks and Corey"
    # Soil Matric potential - "Brooks and Corey"
        alpha = 2 # m-1
        n = 5
        m = 1 - 1/n 
    end

    # How much water
    theta_water = state.θ + state.θi
# @show theta_water

    # Get augmented liquid
    theta_l = augmented_liquid(porosity,S_s,aux.ψ,theta_water) 

    # Get effective saturation
    S_l = effective_saturation(porosity,theta_l)   # 0.2

    # Get matric potential
    ψ_m = matric_potential(flag,alpha,S_l,n,m)

    # This function calculates pressure head ψ of a soil
    aux.ψ = pressure_head(ψ_m,S_l,porosity,S_s,theta_l)  
  
    # Get hydraulic head    
    aux.h = hydraulic_head(aux.z,aux.ψ)        

end

# ---------------- 4b) Calculate state and derivative of T

# Calculate T based on internal energy state variable
function gradvariables!(
    m::SoilModelFull,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)

# HEAT
  #aux.T = state.ρcT / m.ρc(state, aux, t) # TODO: figure out why can't use aux.T here
    # aux.T = temperature_calculator(m.ρc(state, aux, t),state.ρcT,state.θi)
  transform.T = temperature_calculator(m.ρc(state, aux, t),state.ρcT,state.θi)

# WATER
   # Soil Matric potential - "van Genuchten"
    if flag == "van Genuchten"
        alpha = 2 # m-1
        n = 5
        m = 1 - 1/n 
    elseif flag == "Brooks and Corey"
    # Soil Matric potential - "Brooks and Corey"
        alpha = 2 # m-1
        n = 5
        m = 1 - 1/n 
    end

    # How much water
    theta_water = state.θ + state.θi
# @show theta_water

    # Get augmented liquid
    theta_l = augmented_liquid(porosity,S_s,aux.ψ,theta_water) 

    # Get effective saturation
    S_l = effective_saturation(porosity,theta_l)   # 0.2

    # Get matric potential
    ψ_m = matric_potential(flag,alpha,S_l,n,m)

    # This function calculates pressure head ψ of a soil
    aux.ψ = pressure_head(ψ_m,S_l,porosity,S_s,theta_l)  
  
    # Get hydraulic head    
    transform.h = hydraulic_head(aux.z,aux.ψ)    
end

# Gradient of T calculation
function diffusive!(
    m::SoilModelFull,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
  # HEAT
  diffusive.∇T = ∇transform.T
  # WATER
  diffusive.∇h = ∇transform.h
end
# Calculate thermal flux (non-diffusive (?))
function flux_nondiffusive!(
    m::SoilModelFull,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
end
# Calculate thermal flux (diffusive (?))
function flux_diffusive!(
    m::SoilModelFull,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
  )
  # HEAT
  flux.ρcT -= m.κ(state, aux, t) * diffusive.∇T
  # WATER
  flux.θ -= m.K_s(state, aux, t) * diffusive.∇h
   if aux.z == 0
    #@show   aux.T flux.ρcT
    end
end

# ---------------- 4c) Extra Sources

# Introduce sources of energy (e.g. Metabolic heat from microbes)
function source!(
    m::SoilModelFull,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)

# @show state.θi
# @show state.θ
# @show soil_T

# Update sources for ice and liquid
if state.θi > 0
  # Convert liquid water to ice by freezing (or vice versa)
  F_T = calculate_frozen_water(state.θ,state.θi,soil_T)
  # @show F_T
else 
  F_T = 0;
  # @show F_T
end
  
# Source of ice
source.θi = F_T/917 # rho_i = 0.917 # g cm-3, density of ice
source.θ = -F_T/997 # rho_l = 0.997 # g cm-3, density of water

# @show source.θi
# @show source.θ

end

# ---------------- 4d) Initialization

# Initialize z-Profile
function init_aux!(m::SoilModelFull, aux::Vars, geom::LocalGeometry)
  aux.z = geom.coord[3]
  aux.T = m.initialT(aux, 0)
  aux.h = m.initialh(aux, 0) #aux.z+m.initialθ(state, aux, t) #^(-1/0.378))*(-0.3020)
  aux.ψ = m.initialψ(aux, 0)
end
# Initialize State variables from T to internal energy
function init_state!(m::SoilModelFull, state::Vars, aux::Vars, coords, t::Real)
  state.ρcT = m.ρc(state, aux, t) * aux.T
  state.θ = m.initialθ(aux, 0)
  state.θi = m.initialθi(aux, 0)
end

# ---------------- 4e) Boundary Conditions

# Boundary condition function
function boundary_state!(nf, m::SoilModelFull, state⁺::Vars, aux⁺::Vars,
                         nM, state⁻::Vars, aux⁻::Vars, bctype, t, _...)
  if bctype == 1
    # surface temp
    state⁺.ρcT = m.ρc(state⁻, aux⁻, t) * m.surfaceT(state⁻, aux⁻, t)
    # surface water
    state⁺.θ= m.surfaceθ(state⁻, aux⁻, t)
  elseif bctype == 2
    # bottom
    nothing
  end
end
# Boundary condition function - repeated?
function boundary_state!(nf, m::SoilModelFull, state⁺::Vars, diff⁺::Vars,
                         aux⁺::Vars, nM, state⁻::Vars, diff⁻::Vars, aux⁻::Vars,
                         bctype, t, _...)
  if bctype == 1
    # surface temp
    state⁺.ρcT = m.ρc(state⁻, aux⁻, t) * m.surfaceT(state⁻, aux⁻, t)
    # surface water
    state⁺.θ = m.surfaceθ(state⁻, aux⁻, t)
  elseif bctype == 2
    # bottom
    diff⁺.∇T = -diff⁻.∇T
  end
end