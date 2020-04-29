function energy_conservation(x,y,z)
# ------------------------------------------------------
# Input
#   mineral_properties       ! 'Sand','Clay'
#   theta_liq                ! fraction of water that is liquid
#   theta_ice                ! fraction of water that is ice
# ------------------------------------------------------
# Output
#   error                    ! Error
# ------------------------------------------------------
function indefinite_stack_integral!(
    dg::DGModel,
    m::BalanceLaw,
    Q::MPIStateArray,
    auxstate::MPIStateArray,
    t::Real,
    elems::UnitRange = dg.grid.topology.elems,
)
 #vars_integrals,
 #vars_reverse_integrals,
 #indefinite_stack_integral!,
 #reverse_indefinite_stack_integral!,
 #integral_load_aux!,
 #integral_set_aux!


 # Change in energy over one time step  in entire column, from previous time step to following time step (looping over number of layers)
 # cv is heat capactity, dz is thickness of layer i, tsoi is temperature of soil at layer i,
 # tsoi0 is initial temperature of soil at layer i at current time step
 edif = edif + soilvar.cv(i) * soilvar.dz(i) * (soilvar.tsoi(i) - tsoi0(i)) / dt

 # Flux (energy coming into top layer) over one time step
 # thermal conductivity is tk, tsurf is temperature at surface given by boundary condition,
 # tsoi(1) is temperature just below surface layer layer, z is depth of first layer
 soilvar.gsoi = soilvar.tk(1) * (tsurf - soilvar.tsoi(1)) / (0 - soilvar.z(1))

# What I need to know is :
# - where temperature values (equiv of tsoi are being kept), I need previous and current T profiles
# - where initial temperature values are kept
# - where flux is stored (energy coming in from boundary conditions)

# integrate over entire temperature profile at tsoi0
# integrate over entire temperature profile at end of run
# integrate over all energy fluxes

 err = edif - soilvar.gsoi
 if (abs(err) > 1e-03)
  error ('Soil temperature energy conservation error')
 end
 return err
end


