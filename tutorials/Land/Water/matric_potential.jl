# matric_potential.jl: This function calculates matric potential ψ_m of a soil
function matric_potential(flag,alpha,S_l,n,m)

# ------------------------------------------------------
# Input
#   flag                     ! 'van Genuchten','Brooks and Corey'
#   alpha                    ! van Genuchten - inverse reference potential (m-1)
#   S_l                      ! augmented liquid fraction
#   n                        ! van Genuchten - exponent
#   m = 1-1/n                ! van Genuchten - exponent
# ------------------------------------------------------
# Output
#   ψ_m                      ! Soil matric potential
# ------------------------------------------------------

    # Formula for matric potential: ψ_m(S_l) = -alpha^-1 * S_l^(-1/(n*m)) * [1-S_l^(1/m)]^(1/n)

    # Theta(T) = temperature dependence on conductivity
    if flag == "van Genuchten"
        ψ_m = -alpha^-1 * S_l^(-1/(n*m)) * (1-S_l^(1/m))^(1/n)
    elseif flag == "Brooks and Corey"
        ψ_m = -alpha * S_l^(-m)
    else
        ψ_m = 0
    end

    return ψ_m
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov