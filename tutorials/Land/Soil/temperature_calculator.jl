# temperature_calculator.jl: This function calculates soil temperature
function temperature_calculator(c_s,I,theta_ice) 

# ------------------------------------------------------
# Input
#   c_s                      ! Soil heat capacity
#   I                        ! Internal energy of soil
#   theta_ice                ! fraction of water that is ice
# ------------------------------------------------------
# Output
#   T_soil             ! Change in Soil temperature
# ------------------------------------------------------   
    
    # Formula for the temperature soil (K):
    # T = T0 + [ ( I + theta_ice*density_ice*Lf_0 ) / c_s ]

    # Freezing point of water
    T0 = 273.15 # K
    
    # Specific Latent heat of fusion
    Lf_0 = 333.6e3 # J kg-1
    
    # Density of ice
    density_ice = 917 # kg m-3    
        
    # Temperature of soil (K):
    T_soil = T0 + ( ( I + theta_ice*density_ice*Lf_0 ) / c_s )
    
    return T_soil
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov