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

    # How much water
    theta_water = theta_liq + theta_ice

    # Theta(T) = temperature dependence on conductivity
    Theta_T = temperature_dependence(soil_T, soil_Tref)

    # Gamma(θ_i) = frozen soil imdepence factor
    Gamma_thetai = frozen_impedence_factor(theta_ice, porosity)

    # K_sat: Currently using values from Bonan's sp_08_01.m; Need to switch to Global tabulated values from Dai et al. (2019a) [ Sand: K_sat = 10e-2 m s-1 ; Clay: K_sat = 10e-7 m s-1 ]
    if mineral_properties == "Sand"
        K_sat  = 34 / (3600*100)
    elseif mineral_properties == "Clay"
        K_sat  = 0.0443 / (3600*100)
    else
        K_sat  = 1e-3
    end

    # Get augmented liquid
    theta_l = augmented_liquid(porosity,S_s,ψ,theta_water)

    # Get effective saturation
    S_l = effective_saturation(porosity,theta_l)

    # Soil Matric potential - "van Genuchten": Currently using values from Bonan's sp_08_01.m;
    if flag == "van Genuchten"
        alpha = 2.6 # m-1
        n = 1.43
        m = 1 - 1/n
    elseif flag == "Brooks and Corey"
    # Soil Matric potential - "Brooks and Corey"
        alpha = 2.6 # m-1
        n = 1.43
        m = 1 - 1/n
    end

    # Hydraulic conductivity of soil
    if flag == "van Genuchten"
        K_s = Theta_T * Gamma_thetai * K_sat*( S_l^(1/2) * ( 1 - ( 1 - S_l^(1/m) )^m )^2 )
    elseif flag == "Brooks and Corey"
        if S_l < 1
            M = 1/m-1
            K_s = Theta_T * Gamma_thetai * K_sat*( S_l^(2*M+3) )
        elseif S_l >= 1
            K_s = Theta_T * Gamma_thetai * K_sat
        end
    end

    return K_s
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov