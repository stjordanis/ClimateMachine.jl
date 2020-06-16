vars_integrals(::SoilModelMoisture,FT) =  @vars(a::FT)

function vars_state_auxiliary(::SoilModelMoisture, FT)
    @vars begin
        z::FT
        h::FT
        κ::FT
        int::vars_integrals(m, FT)
        a::FT
    end
end


#function update_auxiliary_state!(
#    dg::DGModel,
#    m::SoilModelMoisture,
#    Q::MPIStateArray,
#    t::Real,
#    elems::UnitRange,
#)
#    nodal_update_auxiliary_state!(soil_nodal_update_aux!, dg, m, Q, t, elems)
#    indefinite_stack_integral!(dg, m, Q, dg.state_auxiliary, t, elems)
#  return true
#end

#other stuff we need
function integral_load_auxiliary_state!(
    m::SoilModelMoisture,
    integrand::Vars,
    state::Vars,
    aux::Vars,
)
    integrand.a = state.ν
end

function integral_set_auxiliary_state!(
    m::SoilModelMoisture,
    aux::Vars,
    integral::Vars,
)
    aux.int.a = integral.a
end
    
