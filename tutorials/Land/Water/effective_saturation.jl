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

    return S_l
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov
