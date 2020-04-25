# soil_water_properties.jl: This function calculates soil hydraulic conductivity
function soil_water_properties(mineral_properties,soil_T,soil_Tref,theta_liq,theta_ice,porosity,ψ,S_s,flag) 

# ------------------------------------------------------
# Input
#   mineral_properties       ! 'Sand','Clay'
#   soil_T                   ! Soil temperature
#   soil_Tref                ! Soil reference temperature: annual mean temperature of site
#   theta_liq                ! fraction of water that is liquid
#   theta_ice                ! fraction of water that is ice
#   porosity                 ! Porosity of soil
#   ψ                        ! Soil pressure head 
#   S_s                      ! aquifer specific storage
# ------------------------------------------------------
# Output
#   K_s                      ! Soil hydraulic conductivity
# ------------------------------------------------------   
    
    # Formula for hydraulic conductivity: K_s = Theta(T) * Gamma( θ_i ) * K_sat*{ S_l^1/2 * [ 1 - ( 1 - S_l^(1/m) )^m ]^2 }
            
    # Theta(T) = temperature dependence on conductivity
    Theta_T = temperature_dependence(soil_T, soil_Tref)     
    
    # Gamma(θ_i) = frozen soil imdepence factor
    Gamma_thetai = frozen_impedence_factor(theta_ice, porosity)     
               
    # K_sat: Global tabulated values from Dai et al. (2019a)
    # [ Sand: K_sat = 10e-2 m s-1 ; Clay: K_sat = 10e-7 m s-1 ]
    if mineral_properties == "Sand"
        K_sat  = 1e-3
    elseif mineral_properties == "Clay"
        K_sat  = 1e-7
    else
        K_sat  = 1e-4
    end
    
    # Get augmented liquid
    theta_l = augmented_liquid(porosity,S_s,ψ,theta_liq) 
    
    # Get effective saturation
    S_l = effective_saturation(porosity,theta_l)  

    # Soil Matric potential - "van Genuchten"
    if flag == "van Genuchten"
        alpha = 2 # m-1
        n = 5
        m = 1 - 1/n 
    elseif flag == "Brooks and Corey"
    # Soil Matric potential - "Brooks and Corey"
        alpha = 2 # m-1
        n = 5
        m = 1 - 1/n 
    end
    
    # Matric potential ψ_m of soil    
    ψ_m = matric_potential(flag,alpha,S_l,n,m)
    
    # Pressure head ψ of soil
    ψ = pressure_head(ψ_m,S_l,porosity,S_s,theta_l)  
    
    # Hydraulic conductivity of soil    
    if flag == "van Genuchten"
        K_s = Theta_T * Gamma_thetai * K_sat*( S_l^(1/2) * ( 1 - ( 1 - S_l^(1/m) )^m )^2 )  
    elseif flag == "Brooks and Corey"
        if S_l < 1
            K_s = Theta_T * Gamma_thetai * K_sat*( S_l^(2*m+3) )  
        elseif S_l >= 1
            K_s = Theta_T * Gamma_thetai * K_sat * (1) 
        end
    end
    
    return K_s
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov