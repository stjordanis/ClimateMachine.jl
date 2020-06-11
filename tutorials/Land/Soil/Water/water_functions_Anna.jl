abstract type AbstractWater end
Base.@kwdef struct waterfunctions{HC, MP} <: AbstractWater
    hydraulic_cond::HC     = vanGenuchten{FT}()
    matric_pot::MP = vanGenucthen{FT}()
end

# Three possible models to use for matric potential and hydraulic conductivity.
# Specify parameters of them here.

abstract type AbstractHydraulicsModel end
#Base.@kwdef struct vanGenuchten{FT} <: AbstractHydraulicsModel
#    "n - exponent; that then determines the exponent m used in the model."
#    n::FT = FT(1.43);
#    m::FT = FT(1.0-1.0/n);
#    " alpha  - inverse of this carries units in the expression for matric potential "
#    α::FT = FT(2.6) # inverse meterse
#end

struct vanGenuchten{FT} <: AbstractHydraulicsModel
    "Yolo light clay"
    "n - exponent; that then determines the exponent m used in the model."
     " alpha  - inverse of this carries units in the expression for matric potential (specify in inverse meters)"
    n::FT 
    α::FT
    m::FT
    function vanGenuchten{FT}(;n::FT = FT(1.43), α::FT = FT(2.6))
        new(n, α, FT(1.0-1.0/n))
    end
end

Base.@kwdef struct BrooksCorey{FT} <: AbstractHydraulicsModel
    "m - exponent; ψb - units. Slightly fudged m to better match Havercamp and VG at α=2.0 and n = 2.1."
    ψb::FT = FT(0.1656);# in meters
    m::FT = FT(0.5); #
end

Base.@kwdef struct Haverkamp{FT} <: AbstractHydraulicsModel
    "Yolo light clay"
    "exponent"
    k::FT = FT(1.77);
    "constant A"
    A::FT = FT(124.6/100.0^k) # cm^k. Our sim is in meters - convert
    "constant B"
    B::FT = FT(124.6/100.0^k) # cm^k. Our sim is in meters - convert.
end

##General functions

function hydraulic_head(z,ψ)
    h = z + ψ
    return h
end

function pressure_head(mod::AbstractHydraulicsModel, S_l,porosity,S_s,theta_l)
    if S_l < 1
        ψ = matric_potential(mod, S_l)
    else
        ψ = (theta_l - porosity) / S_s
    end
    return ψ
end

function effective_saturation(porosity,theta_l)
    S_l = theta_l / porosity
    return S_l
end

" 
van Genuchten expression for hydraulic conductivity

"
function hydraulic_conductivity(mod::vanGenuchten, K_sat::FT, S_l::FT, head::FT, z::FT)
  #  S_l(ν) = effective saturation
  #  K_sat = carries the units, constant
  #  head = hydraulic head, units of length. Also a function of state, aux
    #  z  = vertical coordinate, units of length - function of aux
    @unpack n, m = mod;
    FT = typeof(K_sat)
    if S_l < 1.0
        if S_l <= 0.0
            K = FT(0.0)
        else
            K = K_sat*sqrt(S_l)*(1.0-(1.0-S_l^(1.0/m))^m)^2.0
        end
    else
        K = K_sat
    end
    return K
end


"
Brooks and Corey expression for hydraulic conductivity

"
function hydraulic_conductivity(mod::BrooksCorey, K_sat::FT, S_l::FT, head::FT, z::FT)
  #  S_l(ν) = effective saturation
  #  K_sat = carries the units, constant
  #  head = hydraulic head, units of length. Also a function of state, aux
    #  z  = vertical coordinate, units of length - function of aux
    @unpack ψb, m = mod
    FT = typeof(K_sat)
    if S_l < 1.0
        if S_l <= 0.0
            K = FT(0.0)
        else
            K = K_sat*S_l^(2.0*mod.m+3.0)
        end
    else
        K = K_sat
    end
    return K
end

" 
Haverkamp expression for hydraulic conductivity

"
function hydraulic_conductivity(mod::Haverkamp, K_sat::FT, S_l::FT, head::FT, z::FT)
  #  S_l(ν) = effective saturation
  #  K_sat = carries the units, constant
  #  head = hydraulic head, units of meters. Also a function of state, aux
    #  z  = vertical coordinate, units of meters - function of aux
    @unpack k, A, B = mod
    FT = typeof(K_sat)
    if S_l<1
        if S_l <= 0
            K = FT(0.0)
        else
            K = K_sat*A/(B+abs(head-z)^k)
        end
    else
        K = K_sat
    end
    return K
end

" 
van Genuchten expression for matric potential. Also to be used with Haverkamp conductivity

"
function matric_potential(mod::vanGenuchten, S_l::FT)
  #  S_l(ν) = effective saturation
    @unpack n, m, α = mod;
    FT = typeof(S_l)
    if S_l <=0
        ψ_m = -1e30;
    elseif S_l <= 1
        ψ_m = -((S_l^(-1 / m)-1) * α^(-n))^(1 / n)
        #ψ_m = (-alpha^-1 * S_l^(-1/(n*M)) * (1-S_l^(1/M))^(1/n))
    else
        ψ_m = 0
    end
    return ψ_m
end

" 
Brooks and Corey expression for matric potential

"
function matric_potential(mod::BrooksCorey, S_l::FT)
  #  S_l(ν) = effective saturation. This needs to be confirmed.
    @unpack ψb, m = mod;
    FT = typeof(S_l)
    if S_l <=0
        ψ_m = -1e30;
    elseif S_l <= 1
        ψ_m = -ψb*S_l^(-1.0/m)
    else
        ψ_m = ψb
    end
    return ψ_m
end

# "
# Conversion of liquid water to ice by freezing

# "
# function calculate_frozen_water(theta_liq,theta_ice,T)
#     tao_FT = 2e2 #  = max( dt , CFL_bound ), dt = 100 (s) ,  CFL_bound = 200 (s)
#     F_T = ( rho_l*theta_liq*heaviside( (T_f-T) )  - rho_i*theta_ice*heaviside( (T-T_f) ) ) / tao_FT
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
