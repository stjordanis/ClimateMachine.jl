#### Water function

abstract type AbstractWater{FT<:AbstractFloat} end
abstract type AbstractHydraulicsModel{FT<:AbstractFloat} end

const AHM = AbstractHydraulicsModel

"""
    waterfunctions{FT, HC, MP} <: AbstractWater{FT}

Please document
"""
Base.@kwdef struct waterfunctions{FT, HC, MP} <: AbstractWater{FT}
    hydraulic_cond::HC = vanGenuchten{FT}()
    matric_pot::MP = vanGenuchten{FT}()
end

function waterfunctions(;
    hydraulic_cond::AHM{FT} = vanGenuchten{FT}(),
    matric_pot::AHM{FT} = vanGenuchten{FT}()
    ) where {FT}
    return waterfunctions{FT,
        typeof(hydraulic_cond),
        typeof(matric_pot)
        }(hydraulic_cond, matric_pot)
end


# Three possible models to use for matric potential and hydraulic conductivity.
# Specify parameters of them here.

#Base.@kwdef struct vanGenuchten{FT} <: AbstractHydraulicsModel
#    "n - exponent; that then determines the exponent m used in the model."
#    n::FT = FT(1.43);
#    m::FT = FT(1.0-1.0/n);
#    " alpha  - inverse of this carries units in the expression for matric potential "
#    α::FT = FT(2.6) # inverse meterse
#end

"""
    vanGenuchten{FT} <: AbstractHydraulicsModel{FT}

parameters for Yolo light clay
"""
struct vanGenuchten{FT} <: AbstractHydraulicsModel{FT}
    "n - exponent; that then determines the exponent m used in the model."
    n::FT
    "inverse of this carries units in the expression for matric potential (specify in inverse meters)"
    α::FT
    "Exponent parameter"
    m::FT
    function vanGenuchten{FT}(;n::FT = FT(1.43), α::FT = FT(2.6))
        new(n, α, 1-1/FT(n))
    end
end

"""
    BrooksCorey{FT} <: AbstractHydraulicsModel{FT}

Please document

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct BrooksCorey{FT} <: AbstractHydraulicsModel{FT}
    "m - exponent; ψb - units. Slightly fudged m to better match Havercamp and VG at α=2.0 and n = 2.1. (meters)"
    ψb::FT = FT(0.1656);
    "Please document"
    m::FT = FT(0.5);
end

"""
    Haverkamp{FT} <: AbstractHydraulicsModel{FT}

Please document

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Haverkamp{FT} <: AbstractHydraulicsModel{FT}
    "exponent for Yolo light clay"
    k::FT = FT(1.77);
    "constant A cm^k. Our sim is in meters - convert"
    A::FT = FT(124.6/100.0^k)
end

##General functions

hydraulic_head(z,ψ) = z + ψ

"""
    pressure_head(
            model::AbstractHydraulicsModel{FT},
            S_l::FT,
            porosity::FT,
            S_s::FT,
            θ_l::FT
        ) where {FT}

Please document
"""
function pressure_head(
        model::AbstractHydraulicsModel{FT},
        S_l::FT,
        porosity::FT,
        S_s::FT,
        θ_l::FT
    ) where {FT}
    if S_l < 1
        ψ = matric_potential(model, S_l)
    else
        ψ = (θ_l - porosity) / S_s
    end
    return ψ
end

effective_saturation(porosity::FT, θ_l::FT) where {FT} = θ_l / porosity

"
    hydraulic_conductivity(
            model::vanGenuchten{FT},
            K_sat::FT,
            S_l::FT,
            ψ::FT
        ) where {FT}

van Genuchten expression for hydraulic conductivity
"
function hydraulic_conductivity(
        model::vanGenuchten{FT},
        K_sat::FT,
        S_l::FT,
        ψ::FT
    ) where {FT}
  #  S_l(ν) = effective saturation
  #  K_sat = carries the units, constant
  #  ψ = pressure head, units of length. Also a function of state, aux
    @unpack n, m = model;

    if S_l < 1
        if S_l <= 0
            K = FT(0)
        else
            K = K_sat*sqrt(S_l)*(1-(1-S_l^(1/m))^m)^2
        end
    else
        K = K_sat
    end
    return K
end


"
    hydraulic_conductivity(
            model::BrooksCorey{FT},
            K_sat::FT,
            S_l::FT,
            ψ::FT
        ) where {FT}

Brooks and Corey expression for hydraulic conductivity
"
function hydraulic_conductivity(
        model::BrooksCorey{FT},
        K_sat::FT,
        S_l::FT,
        ψ::FT#
        ) where {FT}
  #  S_l(ν) = effective saturation
    #  K_sat = carries the units, constant
    #  ψ = pressure head, units of length. Also a function of state, aux
    @unpack ψb, m = model

    if S_l < 1
        if S_l <= 0
            K = FT(0)
        else
            K = K_sat*S_l^(2 * m + 3)
        end
    else
        K = K_sat
    end
    return K
end

"
    hydraulic_conductivity(
            model::Haverkamp{FT},
            K_sat::FT,
            S_l::FT,
            ψ::FT,
        ) where {FT}

Haverkamp expression for hydraulic conductivity
"
function hydraulic_conductivity(
        model::Haverkamp{FT},
        K_sat::FT,
        S_l::FT,
        ψ::FT,
    ) where {FT}

  #  S_l(ν) = effective saturation
    #  K_sat = carries the units, constant
        #  ψ = pressure head, units of length. Also a function of state, aux
    @unpack k, A = model

    if S_l<1
        if S_l <= 0
            K = FT(0)
        else
            #K = K_sat*A/(B+abs(head-z)^k)
            K = K_sat*A/(A+abs(ψ)^k)
        end
    else
        K = K_sat
    end
    return K
end

"
    matric_potential(
            model::vanGenuchten{FT},
            S_l::FT
        ) where {FT}

van Genuchten expression for matric potential. Also to be used with Haverkamp conductivity
"
function matric_potential(
        model::vanGenuchten{FT},
        S_l::FT
    ) where {FT}
  #  S_l(ν) = effective saturation
    @unpack n, m, α = model;

    S_l = max(S_l, FT(0))
    if S_l <=0
        ψ_m = -1e30;
    elseif  S_l <= 1
        ψ_m = -((S_l^(-1 / m)-1) * α^(-n))^(1 / n)
        #ψ_m = (-alpha^-1 * S_l^(-1/(n*M)) * (1-S_l^(1/M))^(1/n))
    else
        ψ_m = 0
    end
    return ψ_m
end

"
    matric_potential(
            model::BrooksCorey{FT},
            S_l::FT
        ) where {FT}

Brooks and Corey expression for matric potential
"
function matric_potential(
        model::BrooksCorey{FT},
        S_l::FT
    ) where {FT}
    #  S_l(ν) = effective saturation. This needs to be confirmed.
    @unpack ψb, m = model;

    if S_l <=0
        ψ_m = -1e30;
    elseif S_l <= 1
        ψ_m = -ψb*S_l^(-1/m)
    else
        ψ_m = ψb
    end
    return ψ_m
end

# "
# Conversion of liquid water to ice by freezing

# "
# function calculate_frozen_water(θ_liq,theta_ice,T)
#     tao_FT = 2e2 #  = max( dt , CFL_bound ), dt = 100 (s) ,  CFL_bound = 200 (s)
#     F_T = ( rho_l*θ_liq*heaviside( (T_f-T) )  - rho_i*theta_ice*heaviside( (T-T_f) ) ) / tao_FT
#     return F_T
# end

# "
# Expression for frozen soil impedence factor

# "
# function frozen_impedence_factor(theta_ice, porosity)
#     S_i = theta_ice / porosity
#     Gamma_thetai = 10^(-Omega*S_i)
#     return Gamma_thetai
# end

# "
# Heaviside function

# "
# function heaviside(t)
#      Hside = 0.5 * (sign(t) + 1)
#     return Hside
# end

# "
# Expression for temperature dependence of hydraulic conductivity

# "
# function temperature_dependence(soil_T, soil_Tref)
#     Theta_T = exp( gamma_temp*( soil_T - soil_Tref) )
#     return Theta_T
# end
