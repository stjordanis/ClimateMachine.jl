" van Genuchten expression for matric potential. Also to be used with Havercamp conductivity. "
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

" Brooks and Corey expression for matric potential:"
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
