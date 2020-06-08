" van Genuchten expression for hydraulic conductivity:"
function hydraulic_conductivity!(mod::vanGenuchten, K_sat::FT, S_l::FT, head::FT, z::FT)
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


" Brooks and Corey expression for hydraulic conductivity:"
function hydraulic_conductivity!(mod::BrooksCorey, K_sat::FT, S_l::FT, head::FT, z::FT)
  #  S_l(ν) = effective saturation
  #  K_sat = carries the units, constant
  #  head = hydraulic head, units of length. Also a function of state, aux
    #  z  = vertical coordinate, units of length - function of aux
    @unpack n,m = mod
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

" Havercamp expression for hydraulic conductivity:"
function hydraulic_conductivity!(mod::Havercamp, K_sat::FT, S_l::FT, head::FT, z::FT)
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
