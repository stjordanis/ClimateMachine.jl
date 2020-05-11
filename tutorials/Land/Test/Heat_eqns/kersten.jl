# kersten.jl: This function calculates kersten number
function kersten(theta_ice,vom,porosity,Sr,a,b,v_sand,v_gravel)

# ------------------------------------------------------
# Input
#   theta_ice                ! fraction of water that is ice
#   vom                      ! Volume fraction of organic matter in soil
#   porosity                 ! Soil Porosity
#   Sr                       ! Relative Saturation Sr = (theta_liquid + theta_ice) / porosity
#   a                        ! a = -0.24 +/- 0.04 ... adjustable parameter based on soil measurements
#   b                        ! b = 18.1 +/- 1.1 ... adjustable parameter based on soil measurements
#   v_sand                   ! fraction sand
#   v_gravel                 ! fraction gravel
# ------------------------------------------------------
# Output
#   Ke         ! Kersten number
# ------------------------------------------------------

    # If frozen
    if theta_ice > 0
        Ke = Sr^(1+vom)
    else # If not frozen
        Ke = Sr^(0.5*(1+vom- a*v_sand - v_gravel))*( (1 + exp(-b*Sr))^(-3) - ((1-Sr)/2)^3 )^(1-vom)
    end

    return Ke
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov