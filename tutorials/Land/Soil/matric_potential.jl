# matric_potential.jl: This function calculates matric potential ψ_m of a soil
function matric_potential(flag,S_l) 

# ------------------------------------------------------
# Input
#   flag                     ! 'van Genuchten','Brooks and Corey'
#   S_l                      ! augmented liquid fraction
# ------------------------------------------------------
# Output
#   ψ_m                      ! Soil matric potential
# ------------------------------------------------------   
    
    # Formula for matric potential: ψ_m(S_l) = -alpha^-1 * S_l^(-1/(n*m)) * [1-S_l^(1/m)]^(1/n)
    
    # Soil Matric potential - "van Genuchten"
    #   n                        ! van Genuchten - exponent
    #   m = 1-1/n                ! van Genuchten - exponent
    #   alpha                    ! van Genuchten - inverse reference potential (m-1)
    if flag == "van Genuchten"
        alpha = 0.02 # m-1
        n = 5
        M = 1 - 1/n 
    elseif flag == "Brooks and Corey"
    # Soil Matric potential - "Brooks and Corey"
        alpha = 0.02 # m-1
        n = 5
        M = 1 - 1/n 
    end
    
    # Theta(T) = temperature dependence on conductivity
    if flag == "van Genuchten"
        ψ_m = -alpha^-1 * S_l^(-1/(n*M)) * (1-S_l^(1/M))^(1/n) 
    elseif flag == "Brooks and Corey"
        M_b = 1/M-1
        ψ_b = alpha^-1
        ψ_m = -ψ_b * S_l^(-M_b)
    end
    
    return ψ_m
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov