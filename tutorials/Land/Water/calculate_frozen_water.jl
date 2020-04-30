# calculate_frozen_water.jl: This function calculates conversion of liquid water to ice by freezing
function calculate_frozen_water(theta_liq,theta_ice,T) 

# ------------------------------------------------------
# Input
#   theta_liq               ! Soil liquid fraction   
#   theta_ice                ! fraction of water that is ice
#   soil_T                   ! Soil temperature
# ------------------------------------------------------
# Output
#   F_T                     ! [kg m-3 s-1] conversion of liquid water to ice by freezing 
# ------------------------------------------------------   
    
    # Update liquid fraction
    rho_l = 997 # g cm-3, density of water
    rho_i = 917 # g cm-3, density of ice
    T_f = 273.15 # K, freezing temperature
    
    tao_FT = 200 # = max( dt , CFL_bound ), dt = 100 (s) ,  CFL_bound = 200 (s)

    # Convert liquid water to ice by freezing
    F_T = ( rho_l*theta_liq*(T_f-T)  - rho_i*theta_ice*(T-T_f) ) / tao_FT

    # Source of ice and water
    sourceθi = F_T/917 # rho_i = 0.917 # g cm-3, density of ice
    sourceθ = -F_T/997 # rho_l = 0.997 # g cm-3, density of water
    
    # If freezing: Cant have more frozen than liquid available, 
    if sourceθi > theta_liq
        F_T = rho_l*theta_liq         
    end
    
    # If melting: Cant have more liquid than frozen available
    if F_T < -theta_ice
        F_T = -rho_i*theta_ice 
    end
 
    return F_T     
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov
