# hydraulic_head.jl: This function calculates hydraulic head  h of a soil
function hydraulic_head(z,ψ) 

# ------------------------------------------------------
# Input
#   z                      ! soil depth
#   ψ                      ! Soil pressure head
# ------------------------------------------------------
# Output
#   h                      ! Soil hydraulic head
# ------------------------------------------------------   
        
    # Soil hydraulic head as function of depth and pressure head
    h = z + ψ
    
    return h 
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov