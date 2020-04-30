# temperature_dependence.jl: This function calculates temperature dependence of hydraulic conductivity
function temperature_dependence(soil_T, soil_Tref) 

# ------------------------------------------------------
# Input
#   soil_T                 ! Soil temperature
#   soil_Tref              ! Soil reference temperature: annual mean temperature of site
# ------------------------------------------------------
# Output
#   Theta_T                ! Temperature dependence of hydraulic conductivity
# ------------------------------------------------------   
        
    # Formula for Temperature dependence: Theta_T = exp [ γ( T - Tref) ]
    
    # Impedence parameter, from Hansson et al. (2004)
    gamma_temp = 2.64e-2 # K-1 , for Tref = 288 K, implies 30% increase in K for a 10K temp increase
    
    # Temperature dependence of hydraulic conductivity 
    Theta_T = exp( gamma_temp*( soil_T - soil_Tref) )
    
    return Theta_T
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov