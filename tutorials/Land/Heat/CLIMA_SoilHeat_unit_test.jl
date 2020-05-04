# --------------------------------- CLIMA SOIL MODEL -----------------------
# CLIMA_SoilHeat.jl: This model simulates soil heat dynamics for the CliMA model

"""
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

# --------------------------------- 1) Import/Export Needed Functions -----------------------

# Add necessary CliMA functions and sub-routines
using StaticArrays
using CLIMA.VariableTemplates
import CLIMA.DGmethods: BalanceLaw,
                        vars_state_auxiliary, vars_state_conservative, vars_state_gradient, vars_state_gradient_flux, vars_integrals, vars_reverse_integrals,
                        flux_first_order!, flux_second_order!, source!,
                        compute_gradient_argument!, compute_gradient_flux!, update_auxiliary_state!, nodal_update_auxiliary_state!, integral_load_auxiliary_state!, integral_set_auxiliary_state!, reverse_integral_load_auxiliary_state!, reverse_integral_set_auxiliary_state!,
                        indefinite_stack_integral!, reverse_indefinite_stack_integral!,
                        init_state_auxiliary!, init_state_conservative!,
                        boundary_state!, wavespeed, LocalGeometry


# --------------------------------- 2) Define Structs ---------------------------------------

# Introduce needed variables into SoilModel struct
Base.@kwdef struct SoilModel{Fρc, Fκ, FiT, Fst} <: BalanceLaw

  # Define heat capacity. This is an input to the model now.
  ρc::Fρc       = (state, aux, t) -> 2.49e6   # [ Sand: ρc = 2.49e6 J m-3 K-1 ; Clay: ρc = 2.61e6 J m-3 K-1 ]
  # Replace this with a function that calculates heat capacity (based on liquid+ice)
  # OR Replace this with tabulated values of heat capacity (based on liquid+ice)

  # Define kappa (thermal conductivity). This is an input to the model now.
  κ::Fκ         = (state, aux, t) -> 2.42     # [ Sand: λ = 2.42 W m-1 K-1 ; Clay: λ = 1.17 W m-1 K-1 ]

  # Define initial and boundary condition parameters
  initialT::FiT = (aux) -> aux.z - (aux.z)^2  # 273.15 # (5*exp(-(aux.z-0.5)^2/(2*(0.2)^2))) # Initial Temperature. This is an input to the model now.
  surfaceT::Fst = (state, aux, t) -> 273.15+10*heaviside(t-t1)-10*heaviside(t-t2) #273.15 # (273.15 + 15.0) + 0.5*10.0 * sinpi(2*(t/(60*60)-8)/24) #(273.15 + 2.0) # Surface boundary condition. This is an input to the model now.
end

# --------------------------------- 3) Define CliMA vars ---------------------------------------

# Stored in the aux state are:
#   `coord` coordinate points (needed for BCs)
#   `u` advection velocity
#   `D` Diffusion tensor
vars_integrals(::SoilModel,FT) = @vars(a::FT) # location to store integrands for bottom up integrals
vars_reverse_integrals(::SoilModel, FT) = @vars(a::FT) # location to store integrands for top down integrals
vars_state_auxiliary(::SoilModel, FT) = @vars(z::FT, T::FT) # stored dg.auxstate
vars_state_conservative(::SoilModel, FT) = @vars(ρcT::FT) #analytical_flux::FT stored in Q , (\rho  c T) is number rows 
#@vars(int::vars_integrals(m, T), rev_int::vars_integrals(m, T), r::T, a::T)
vars_state_gradient(::SoilModel, FT) = @vars(T::FT) # not stored
vars_state_gradient_flux(::SoilModel, FT) = @vars(∇T::SVector{3,FT}) # stored in dg.diffstate

# integrate over entire temperature profile at tsoi0
# integrate over entire temperature profile at end of run
# integrate over all energy fluxes

#add to list of callbacks

# --------------------------------- 4) CliMA functions needed for simulation -------------------
# ---------------- 4a) Update states
# Update all auxiliary variables
function update_auxiliary_state!(
  dg::DGModel,
  m::SoilModel,
  Q::MPIStateArray,
  t::Real,
  elems::UnitRange,
)
  nodal_update_auxiliary_state!(soil_nodal_update_aux!, dg, m, Q, t, elems)
  indefinite_stack_integral!(dg, m, Q, dg.state_auxiliary, t, elems)
  #reverse_indefinite_stack_integral!(dg, m, Q, dg.auxstate, t, elems)
  # Flux (energy coming into top layer) over one time step
  # thermal conductivity is tk, tsurf is temperature at surface given by boundary condition,
  # tsoi(1) is temperature just below surface layer layer, z is depth of first layer
  # soilvar.gsoi = soilvar.tk(1) * (tsurf - soilvar.tsoi(1)) / (0 - soilvar.z(1))
  # analytical_flux=m.ρc(state, aux, t)*m.surfaceT(state, aux, t)
  
  return true
end
#
function integral_load_auxiliary_state!(
    m::SoilModel,
    integrand::Vars,
    state::Vars,
    aux::Vars,
)
    integrand.a = state.ρcT
end

function integral_set_auxiliary_state!(
    m::SoilModel,
    aux::Vars,
    integral::Vars,
)
    aux.int.a = integral.a
end

function reverse_integral_load_auxiliary_state!(
    m::SoilModel,
    integral::Vars,
    state::Vars,
    aux::Vars,
)
    integral.a = aux.int.a
end

function reverse_integral_set_auxiliary_state!(
    m::SoilModel,
    aux::Vars,
    integral::Vars,
)
    aux.rev_int.a = integral.a
end

# Update all auxiliary nodes
function soil_nodal_update_aux!(
  m::SoilModel,
  state::Vars,
  aux::Vars,
  t::Real)
  aux.T = state.ρcT / m.ρc(state, aux, t) # TODO: figure out why can't use aux.T here
end

# ---------------- 4b) Calculate state and derivative of T

# Calculate T based on internal energy state variable
function compute_gradient_argument!(
    m::SoilModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
  transform.T = state.ρcT / m.ρc(state, aux, t)
end
# Gradient of T calculation
function compute_gradient_flux!(
    m::SoilModel,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
  diffusive.∇T = ∇transform.T
end
# Calculate thermal flux (non-diffusive (?))
function flux_first_order!(
    m::SoilModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
end
# Calculate thermal flux (diffusive (?))
function flux_second_order!(
    m::SoilModel,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
  )
   flux.ρcT -= m.κ(state, aux, t) * diffusive.∇T
   if aux.z == 0
    #@show   aux.T flux.ρcT
    end
end

# ---------------- 4c) Extra Sources
# Introduce sources of energy (e.g. Metabolic heat from microbes)
# function source!(m::SoilModel, state::Vars, _...)
# end
function source!(
    m::SoilModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)
#dirac_space=1000000*exp(-(aux.z-(-0.5))^2/(2*(0.1)^2)) #*m.ρc(state, aux, t)/timeend
#dirac_time=1*exp(-(t-(dt))^2/(2*(dt)^2))
#source.ρcT=dirac_space*dirac_time
  # @show(source.ρcT)
end

# ---------------- 4d) Initialization

# Initialize z-Profile
function init_state_auxiliary!(m::SoilModel, aux::Vars, geom::LocalGeometry)
  aux.z = geom.coord[3]
  aux.T = m.initialT(aux)
end
# Initialize State variables from T to internal energy
function init_state_conservative!(m::SoilModel, state::Vars, aux::Vars, coords, t::Real)
  state.ρcT = m.ρc(state, aux, t) * aux.T
end

# ---------------- 4e) Boundary Conditions

# Boundary condition function
function boundary_state!(nf, m::SoilModel, state⁺::Vars, aux⁺::Vars,
                         nM, state⁻::Vars, aux⁻::Vars, bctype, t, _...)
  if bctype == 1
    # surface
    # state⁺.ρcT = m.ρc(state⁻, aux⁻, t) * m.surfaceT(state⁻, aux⁻, t)
    nothing
  elseif bctype == 2
    # bottom
    nothing
  end
end
# Boundary condition function - repeated?
function boundary_state!(nf, m::SoilModel, state⁺::Vars, diff⁺::Vars,
                         aux⁺::Vars, nM, state⁻::Vars, diff⁻::Vars, aux⁻::Vars,
                         bctype, t, _...)
  if bctype == 1
    # surface
    diff⁺.∇T = -diff⁻.∇T
    # state⁺.ρcT = m.ρc(state⁻, aux⁻, t) * m.surfaceT(state⁻, aux⁻, t)
  elseif bctype == 2
    # bottom
    diff⁺.∇T = -diff⁻.∇T
  end
end

