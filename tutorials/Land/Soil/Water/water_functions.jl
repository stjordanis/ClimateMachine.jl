using DocStringExtensions

abstract type AbstractWater{FT<:AbstractFloat} end
abstract type AbstractHydraulicsModel{FT<:AbstractFloat} end

const AHM = AbstractHydraulicsModel

"""
    waterfunctions{FT, HC, MP} <: AbstractWater{FT}

The chosen hydraulic model functions for matric potential and hydraulic conductivity.

# Constructors
    waterfunctions(; hydraulic_cond::AHM{FT} = vanGenuchten{FT}(), matric_pot::AHM{FT} = vanGenuchten{FT}())

Creates instances of the hydraulics models given as keyword arguments. These instances specify the set parameters
of that particular hydraulic model. Options are Haverkamp hydraulic conductivity, van Genuchten conductivity and
matric potential, or the Brooks and Corey hydraulic conductivity and matric potential.

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

"""
    vanGenuchten{FT} <: AbstractHydraulicsModel{FT}

The necessary parameters for the van Genuchten hydraulic model; defaults are for Yolo light clay.

# Fields

$(DocStringExtensions.FIELDS)
"""
struct vanGenuchten{FT} <: AbstractHydraulicsModel{FT}
    "Exponent parameter - using in matric potential"
    n::FT
    "used in matric potential. The inverse of this carries units in the expression for matric potential (specify in inverse meters)."
    α::FT
    "Exponent parameter - determined by n, used in hydraulic conductivity"
    m::FT
    function vanGenuchten{FT}(;n::FT = FT(1.43), α::FT = FT(2.6)) where {FT}
        new(n, α, 1-1/FT(n))
    end
end

"""
    BrooksCorey{FT} <: AbstractHydraulicsModel{FT}

The necessary parameters for the Brooks and Corey hydraulic model.

Defaults are chosen to somewhat mirror the Havercamp/vG Yolo light clay hydraulic conductivity/matric potential.

# Fields

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct BrooksCorey{FT} <: AbstractHydraulicsModel{FT}
    "ψ_b - used in matric potential. Units of meters."
    ψb::FT = FT(0.1656);
    "Exponent used in matric potential and hydraulic conductivity."
    m::FT = FT(0.5);
end

"""
    Haverkamp{FT} <: AbstractHydraulicsModel{FT}

The necessary parameters for the Haverkamp hydraulic model for Yolo light clay.

Note that this only is used in creating a hydraulic conductivity function, and another formulation for matric potential must be used.

# Fields

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Haverkamp{FT} <: AbstractHydraulicsModel{FT}
    "exponent"
    k::FT = FT(1.77);
    "constant A (units of cm^k). Our sim is in meters - convert to meters with factor of 1/100^k."
    A::FT = FT(124.6/100.0^k)
end


"""
    hydraulic_head(z,ψ)

Return the hydraulic head.

The hydraulic head is defined as the sum of vertical height and pressure head; meters.
"""
hydraulic_head(z,ψ) = z + ψ

"""
   effective_saturation(porosity::FT, θ_l::FT)

Compute the effective saturation of soil.

θ_l is defined to be zero or positive. If θ_l is negative, hydraulic functions that take it as an argument will return imaginary numbers, resulting in domain errors. However, it is possible that our current solver returns a negative θ_l due to numerical issues. Provide a warning in this case, and correct the value of θ_l so that the integration can proceed. We will remove this once the numerical issues are resolved.
"""
function effective_saturation(
        porosity::FT,
        θ_l::FT
        ) where {FT}

    if θ_l<0
        @show θ_l
        @warn("Augmented liquid fraction is negative - domain error. Artificially setting equal to zero to proceed. ")
        θ_l = 0
    end
    S_l = θ_l / porosity
    return S_l
end

"""
    pressure_head(
            model::AbstractHydraulicsModel{FT},
            porosity::FT,
            S_s::FT,
            θ_l::FT
        ) where {FT}

Determine the pressure head in both saturated and unsaturated soil.
"""
function pressure_head(
        model::AbstractHydraulicsModel{FT},
        porosity::FT,
        S_s::FT,
        θ_l::FT
        ) where {FT}
    
    S_l = effective_saturation(porosity, θ_l)
    if S_l < 1
        ψ = matric_potential(model, S_l)
    else
        ψ = (θ_l - porosity) / S_s
    end
    return ψ
end



"
    hydraulic_conductivity(
            model::vanGenuchten{FT},
            K_sat::FT,
            S_l::FT,
            ψ::FT
        ) where {FT}

Compute the van Genuchten function for hydraulic conductivity.

the van Genuchten and Brooks and Corey expressions for conductivity require the effective saturation as an argument, while the Haverkamp expression require the pressure head. We've created this function to run using multiple dispatch, so it takes as arguments both pressure head and effective saturation.
"
function hydraulic_conductivity(
        model::vanGenuchten{FT},
        K_sat::FT,
        S_l::FT,
        ψ::FT
    ) where {FT}
    @unpack n, m = model;
    if S_l < 1
        K = K_sat*sqrt(S_l)*(1-(1-S_l^(1/m))^m)^2
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

Compute the Brooks and Corey function for hydraulic conductivity.

the van Genuchten and Brooks and Corey expressions for conductivity require the effective saturation as an argument, while the Haverkamp expression require the pressure head. We've created this function to run using multiple dispatch, so it takes as arguments both pressure head and effective saturation.

"
function hydraulic_conductivity(
        model::BrooksCorey{FT},
        K_sat::FT,
        S_l::FT,
        ψ::FT
        ) where {FT}
    @unpack ψb, m = model

    if S_l < 1
        K = K_sat*S_l^(2 * m + 3)
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

Compute the Haverkamp function for hydraulic conductivity.

the van Genuchten and Brooks and Corey expressions for conductivity require the effective saturation as an argument, while the Haverkamp expression require the pressure head. We've created this function to run using multiple dispatch, so it takes as arguments both pressure head and effective saturation.

"
function hydraulic_conductivity(
        model::Haverkamp{FT},
        K_sat::FT,
        S_l::FT,
        ψ::FT,
    ) where {FT}
    @unpack k, A = model

    if S_l<1
        K = K_sat*A/(A+abs(ψ)^k)
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

Compute the van Genuchten function for matric potential.

This is also to be used with the Haverkamp hydraulic conductivity function.
"
function matric_potential(
        model::vanGenuchten{FT},
        S_l::FT
    ) where {FT}
    @unpack n, m, α = model;

    if S_l < 1
        ψ_m = -((S_l^(-1 / m)-1) * α^(-n))^(1 / n)
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

Compute the Brooks and Corey function for matric potential.
"
function matric_potential(
        model::BrooksCorey{FT},
        S_l::FT
    ) where {FT}
    @unpack ψb, m = model;

    if S_l <= 1
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
