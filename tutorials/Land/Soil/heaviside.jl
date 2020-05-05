# heaviside.jl: This function calculates heaviside function
function heaviside(t) 

# ------------------------------------------------------
# Input
#   t               ! is this value positive or negative?  
# ------------------------------------------------------
# Output
#   Hside           ! if 'yes = 1, if 'no' = 0 
# ------------------------------------------------------   
    
    # Heaviside function
     Hside = 0.5 * (sign(t) + 1)
 
    return Hside     
end

# ______________________________________________________________________________________________________________________
# Written by: Elias Massoud, Jet Propulsion Laboratory/California Institute of Technology, elias.massoud@jpl.nasa.gov
