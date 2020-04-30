# effective_saturation.jl: This function calculates effective saturation S_l of a soil
function effective_saturation(porosity,S_s,ψ,theta_liq) 

# ------------------------------------------------------
# Input
#   porosity                ! soil porosity
#   S_s                     ! aquifer specific storage
#   ψ                       ! Soil pressure head
#   theta_liq               ! Soil liquid fraction   
# ------------------------------------------------------
# Output
#   S_l                      ! soil effective saturation
# ------------------------------------------------------   
    
    # Update liquid fraction
    Σ = S_s * ψ
    
    # theta_l = augmented liquid; integrates specific storage S_s from water table at ψ=0 to saturated location of interest 
    if theta_liq <= porosity
        theta_l = theta_liq
    elseif theta_liq > porosity
        theta_l = theta_liq  + Σ 
    end
    
    # Effective saturation
    S_l = theta_l / porosity
    
    return S_l 
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov
