# pressure_head.jl: This function calculates pressure head ψ of a soil
function pressure_head(ψ_m,S_l,porosity,S_s,theta_l)

# ------------------------------------------------------
# Input
#   ψ_m                      ! soil matric potential
#   S_l                      ! augmented liquid fraction
#   porosity                 ! soil porosity
#   S_s                      ! aquifer specific storage
#   theta_l                  ! soil augmented liquid
# ------------------------------------------------------
# Output
#   ψ                        ! Soil pressure head
# ------------------------------------------------------

    # Soil pressure head as function of saturation level
    if S_l < 1
        ψ = ψ_m
    elseif S_l >= 1
        ψ =  ( theta_l - porosity ) / S_s
    else
        ψ = ψ_m
    end

    return ψ
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov