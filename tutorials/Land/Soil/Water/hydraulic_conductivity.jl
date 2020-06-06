# hydraulic_conductivity.jl: This function calculates hydraulic conductivity.
#Note that we have things defined in multiple places - the parameters of the VG model and BC model are here and in matric potential. we should pass them with the flag so that if they change we dont need to update in multiple places.
function lower_quadratic(a, b, c)
  discr = b^2 - 4*a*c
  discr >= 0 ?   (-b - sqrt(discr))/2a : NaN
end # function

function hydraulic_conductivity(K_sat, S_l, head, zed, flag)

    if flag == "van Genuchten"
        
        n = 2
        M = 1 - 1/n
        #if S_l > 0.995
        #    S_l = 0.995#
        #
        #        end
        if S_l <=0
            K = 0.0;
        elseif S_l < 1.1
            k = sqrt(S_l)*(1.0-(1.0-real(S_l^(1.0/M))^M)^2.0)
            K = lower_quadratic(0.9999,-(K_sat+k),K_sat*k);
        else
            K = K_sat
        end

    elseif flag == "Brooks and Corey"
        alpha = 0.02 # m-1
        n = 5
        M = 1 - 1/n
        if S_l <=0
            K = 0.0;
        elseif S_l < 1.0
            K = K_sat*S_l^(2.0*M+3.0)
        else
            K = K_sat
        end
    elseif flag == "Havercamp"
        K = K_sat*124.6/(124.6+abs(100*(head-zed))^1.77)
        #K = K_sat*124.6/(124.6+abs(100*(head-zed))^1.77)    
    end
    return K
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
