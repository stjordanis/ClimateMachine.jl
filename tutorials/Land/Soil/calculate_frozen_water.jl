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
    rho_l = 997 # kg m-3, density of water
    rho_i = 917 # kg m-3, density of ice
    T_f = 273.15 # K, freezing temperature
    
    tao_FT = 2e6 #  = max( dt , CFL_bound ), dt = 100 (s) ,  CFL_bound = 200 (s)

    # Convert liquid water to ice by freezing
    F_T = ( rho_l*theta_liq*heaviside( (T_f-T) )  - rho_i*theta_ice*heaviside( (T-T_f) ) ) / tao_FT         

    return F_T     
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov
