abstract type AbstractWater end
Base.@kwdef struct waterfunctions{HC, MP} <: AbstractWater
    hydraulic_cond::HC     = vanGenuchten{FT}()
    matric_pot::MP = vanGenucthen{FT}()
end

# Three possible models to use for matric potential and hydraulic conductivity.
# Specify parameters of them here.

abstract type AbstractHydraulicsModel end
Base.@kwdef struct vanGenuchten{FT} <: AbstractHydraulicsModel
    "n - exponent; that then determines the exponent m used in the model."
    n::FT = FT(1.43);
    m::FT = FT(1.0-1.0/n);
    " alpha  - inverse of this carries units in the expression for matric potential "
    α::FT = FT(2.6) # inverse meterse
end

Base.@kwdef struct BrooksCorey{FT} <: AbstractHydraulicsModel
    "m - exponent; ψb - units. TBD the representative value"
    ψb::FT = FT(0.0);
    m::FT = FT(1.0-1.0/n);
end

Base.@kwdef struct Havercamp{FT} <: AbstractHydraulicsModel
    "exponent"
    k::FT = FT(1.77);
    "constant A"
    A::FT = FT(124.6/100.0^k) # these carry units of cm^k. Our sim is in meters - convert
    "constant B"
    B::FT = FT(124.6/100.0^k) # cm^k. Our sim is in meters - convert.
end

##General functions

function hydraulic_head(z,ψ)
# ------------------------------------------------------
# Input
#   z                      ! soil depth [in m]
#   ψ                      ! Soil pressure head [m]
# ------------------------------------------------------
# Output
#   h                      ! Soil hydraulic head [m]
# ------------------------------------------------------
    h = z + ψ
    return h
end

function pressure_head(mod::AbstractHydraulicsModel, S_l,porosity,S_s,theta_l)
# ------------------------------------------------------
# Input
#   ψ_m                      ! soil matric potential [m]
#   S_l                      ! augmented liquid fraction
#   porosity                 ! soil porosity
#   S_s                      ! aquifer specific storage
#   theta_l                  ! soil augmented liquid
# ------------------------------------------------------
# Output
#   ψ                        ! Soil pressure head [m]
# ------------------------------------------------------
    if S_l < 1
        ψ = matric_potential!(mod, S_l)
    else
        ψ = (theta_l - porosity) / S_s
    end
    return ψ
end

function effective_saturation(porosity,theta_l)
# ------------------------------------------------------
# Input
#   porosity                ! soil porosity
#   theta_l                 ! soil augmented liquid
# ------------------------------------------------------
# Output
#   S_l                      ! soil effective saturation
# ------------------------------------------------------
    S_l = theta_l / porosity
    return S_l
end
