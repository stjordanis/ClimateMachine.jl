# hydraulic_conductivity.jl: This function calculates hydraulic conductivity.
#Note that we have things defined in multiple places - the parameters of the VG model and BC model are here and in matric potential. we should pass them with the flag so that if they change we dont need to update in multiple places.
using Parameters 
abstract type AbstractHydraulicsModel end

Base.@kwdef struct vanGenuchten{FT} <: AbstractHydraulicsModel
    "comment here"
    n::FT = FT(2)
    "comment here"
    m::FT = FT(1 - 1/n)
    "comment here"
    α::FT = FT(2)
end


function hydraulic_conductivity(mod::vanGenuchten, K_sat, S_l)
    @unpack n,m = mod;
    FT = typeof(K_sat)
    if S_l <=0
        K = FT(0.0);
    elseif S_l <=1
        K = K_sat * sqrt(S_l) * (1 - (1 - (S_l^(1/m))^m))^2
    else
        K = K_sat
    end
end

function matric_potential(mod::vanGenuchten,S_l)
    @unpack n,m,α = mod;
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

function saturation_from_potential(mod::vanGenuchten,ψ)
    @unpack n,m,α = mod;
    (1 + (α * abs(ψ))^n)^(-m);
end

function hydraulic_conductivity2(K_sat, ψ,  n, α)
    m = 1.0 - 1.0/n
    if (ψ <= 0)
        Se = (1 + (α * abs(ψ))^n)^-m
    else
        Se = 1
    end

    if Se <= 1
        K = K_sat * sqrt(Se) * (1 - (1 - Se^(1/m))^m)^2
    else
        K = K_sat
    end
    if Se > 1
        @show K, Se
    elseif Se<0
        @show K, Se
    end
    return K
end
