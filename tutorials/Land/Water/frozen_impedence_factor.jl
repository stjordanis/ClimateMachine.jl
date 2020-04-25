# frozen_impedence_factor.jl: This function calculates frozen soil Impedence factor
function frozen_impedence_factor(theta_ice, porosity) 

# ------------------------------------------------------
# Input
#   theta_ice                ! fraction of water that is ice
#   porosity                 ! Soil porosity
# ------------------------------------------------------
# Output
#   Gamma_thetai             ! Frozen soil Impedence factor
# ------------------------------------------------------   
    
    # Impedence parameter, from Hansson et al. (2004)
    Omega = 7 
    
    # Impedence factor
    S_i = theta_ice / porosity 
    Gamma_thetai = 10^(-Omega*S_i)
    
    return Gamma_thetai
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov
