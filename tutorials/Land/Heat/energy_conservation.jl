# heat_capacity.jl: This function calculates soil heat capacity
function energy_conservation(x,y,z)
# ------------------------------------------------------
# Input
#   mineral_properties       ! 'Sand','Clay'
#   theta_liq                ! fraction of water that is liquid
#   theta_ice                ! fraction of water that is ice
# ------------------------------------------------------
# Output
#   κ_out                    ! Soil thermal conductivity
# ------------------------------------------------------
 #indefinite_stack_integral!
 edif = edif + soilvar.cv(i) * soilvar.dz(i) * (soilvar.tsoi(i) - tsoi0(i)) / dt
 soilvar.gsoi = soilvar.tk(1) * (tsurf - soilvar.tsoi(1)) / (0 - soilvar.z(1))
 err = edif - soilvar.gsoi
 if (abs(err) > 1e-03)
  error ('Soil temperature energy conservation error')
 end
 return err
end
