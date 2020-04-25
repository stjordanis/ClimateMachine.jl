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
                        vars_aux, vars_state, vars_gradient, vars_diffusive,
                        flux_nondiffusive!, flux_diffusive!, source!,
                        gradvariables!, diffusive!, update_aux!, nodal_update_aux!,
                        init_aux!, init_state!,
                        boundary_state!, wavespeed, LocalGeometry


# --------------------------------- 2) Define Structs ---------------------------------------


# Introduce needed variables into SoilModel struct
Base.@kwdef struct SoilModelMoisture{Fκ, Fiθ, Fsθ} <: BalanceLaw
  # Define kappa (hydraulic conductivity)
  K_s::Fκ         = (state, aux, t) -> (1e-3*(0.34/(60*60))*1.175e6/((1.175e6+abs.(aux.h-aux.z)^4.74))) #(0.001/(60*60*24)) [m/s] typical value taken from Land Surface Model CLiMA, table 2.2, =0.1cm/day (0.34*1.175e6/(1.175+abs.(aux.h)^4.74)) 
  # Define initial and boundary condition parameters
  initialθ::Fiθ = (state, aux, t) -> 0.1 # [m3/m3] constant water content in soil, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
  surfaceθ::Fsθ = (state, aux, t) -> 0.267 #267 # [m3/m3] constant flux at surface, from Bonan, Ch.8, fig 8.8 as in Haverkamp et al. 1977, p.287
end


# --------------------------------- 3) Define CliMA vars ---------------------------------------


# Stored in the aux state are:
#   `coord` coordinate points (needed for BCs)
#   `u` advection velocity
#   `D` Diffusion tensor
vars_aux(::SoilModelMoisture, Fθ) = @vars(z::Fθ, h::Fθ , ψ::Fθ) # p::Fθ stored in dg.auxstate
vars_state(::SoilModelMoisture, Fθ) = @vars(θ::Fθ, θi::Fθ) # stored in Q
vars_gradient(::SoilModelMoisture, Fθ) = @vars(h::Fθ) # not stored
vars_diffusive(::SoilModelMoisture, Fθ) = @vars(∇h::SVector{3,Fθ}) # stored in dg.diffstate


# --------------------------------- 4) CliMA functions needed for simulation -------------------

# ---------------- 4a) Update states

# Update all auxiliary variables
function update_aux!(
    dg::DGModel,
    m::SoilModelMoisture,
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
)
  nodal_update_aux!(soil_nodal_update_aux!, dg, m, Q, t, elems)
  return true
end
# Update all auxiliary nodes
function  soil_nodal_update_aux!(
  m::SoilModelMoisture,
  state::Vars,
  aux::Vars,
  t::Real)
    # flag = "van Genuchten" # - "Brooks and Corey"
    
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
    
    # Get augmented liquid
    theta_l = augmented_liquid(porosity,S_s,ψ,theta_liq) 
    
    # Get effective saturation
    S_l = effective_saturation(porosity,theta_l)   # 0.2
    
    # Get matric potential
    ψ_m = matric_potential(flag,alpha,S_l,n,m)

    # This function calculates pressure head ψ of a soil
    aux.ψ = pressure_head(ψ_m,S_l,porosity,S_s,theta_l)  
       
    # Get hydraulic head      
    aux.h = hydraulic_head(aux.z,aux.ψ)   
    #aux.h = aux.z+((-1/2.7)*(state.θ/1.987)^(-1/3.96))*(1-state.θ/1.987)^(1/3.96)
end


# ---------------- 4b) Calculate state and derivative of theta

# Calculate h based on state variable
function gradvariables!(
    m::SoilModelMoisture,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)     
    # flag = "van Genuchten" # - "Brooks and Corey"
    
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

    # Get augmented liquid
    theta_l = augmented_liquid(porosity,S_s,ψ,theta_liq) 
    #@show theta_l
    # Get effective saturation
    S_l = effective_saturation(porosity,theta_l)   # 0.2
    #@show S_l
    # Get matric potential
    ψ_m = matric_potential(flag,alpha,S_l,n,m)
    #@show ψ_m
    # This function calculates pressure head ψ of a soil
    aux.ψ = pressure_head(ψ_m,S_l,porosity,S_s,theta_l)  
    #@show aux.ψ
    # Get hydraulic head    
    transform.h = hydraulic_head(aux.z,aux.ψ)    
    # transform.h = aux.z+((-1/2.7)*(state.θ/1.987)^(-1/3.96))*(1-state.θ/1.987)^(1/3.96)
end

# Gradient of h calculation
function diffusive!(
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
function  flux_nondiffusive!(
    m::SoilModelMoisture,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
end

# Calculate water flux (diffusive)
function flux_diffusive!(
    m::SoilModelMoisture,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
  )
   flux.θ -= m.K_s(state, aux, t) * diffusive.∇h
    # @show flux.θ
   if aux.z == 0
   end
end

# ---------------- 4c) Extra Sources

# Introduce sources of energy (e.g. Metabolic heat from microbes)
function source!(m::SoilModelMoisture, state::Vars, _...)
end

# ---------------- 4d) Initialization 


# Initialize z-Profile ### what role does this play? when?
function init_aux!(m::SoilModelMoisture, aux::Vars, geom::LocalGeometry)
  aux.z = geom.coord[3]
  aux.h = -100 #aux.z+m.initialθ(state, aux, t) #^(-1/0.378))*(-0.3020)
end

# Initialize State variables from T to internal energy
function init_state!(m::SoilModelMoisture, state::Vars, aux::Vars, coords, t::Real)
  state.θ = m.initialθ(state, aux, t)
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
    nothing
    #diff⁺.∇h = -diff⁻.∇h
    #diff⁺.∇θ = -diff⁻.∇θ        
  end
end
