# hydraulic_conductivity.jl: This function calculates hydraulic conductivity.
#Note that we have things defined in multiple places - the parameters of the VG model and BC model are here and in matric potential. we should pass them with the flag so that if they change we dont need to update in multiple places.
function hydraulic_conductivity(K_sat, S_l, head, zed, flag)

# ------------------------------------------------------
# Input
#   K_sat
#   S_l(ν) = effective saturation
#   hydraulic head
#   z
#   flag                     ! 'van Genuchten','Brooks and Corey', 'Havercamp'

# ------------------------------------------------------
# Output
#   κ                      ! Soil hydraulic conductivity
# ------------------------------------------------------

    # "van Genuchten"
    #   n                        ! van Genuchten - exponent
    #   m = 1-1/n                ! van Genuchten - exponent
    #   alpha                    ! van Genuchten - inverse reference potential (m-1)
    if flag == "van Genuchten"
        n = 1.43
        M = 1.0 - 1.0/n
        if S_l < 1.0
            K = K_sat*sqrt(S_l)*(1.0-(1.0-S_l^(1.0/M))^M)^2.0
        else
            K = 1.0
        end

    elseif flag == "Brooks and Corey"
        alpha = 0.02 # m-1
        n = 5
        M = 1 - 1/n
        if S_1 < 1.0
            K = K_sat*S_l^(2.0*M+3.0)
        else
            K = 1.0
        end
    elseif flag == "Havercamp"
        K = K_sat*124.6/(124.6+abs(100*(head-zed))^1.77)    
    end
    return K
end
