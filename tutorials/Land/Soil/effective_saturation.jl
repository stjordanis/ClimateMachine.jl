# effective_saturation.jl: This function calculates effective saturation S_l of a soil
function effective_saturation(porosity,theta_l) 

# ------------------------------------------------------
# Input
#   porosity                ! soil porosity
#   theta_l                 ! soil augmented liquid   
# ------------------------------------------------------
# Output
#   S_l                      ! soil effective saturation
# ------------------------------------------------------   
    
    # Effective saturation
    S_l = theta_l / porosity
    
#     # Eliminiate super saturation (for now)
#     if S_l >1
#         S_l = 0.9999;
#     end
#     # Eliminiate super dryness (for now)
#     if S_l <0
#         S_l = 0.0001;
#     end
    
    return S_l 
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov
