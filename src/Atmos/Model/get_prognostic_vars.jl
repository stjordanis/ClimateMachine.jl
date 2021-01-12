##### Get prognostic variable list

import ..BalanceLaws: prognostic_vars, get_prog_state

prognostic_vars(::DryModel) = ()
prognostic_vars(::EquilMoist) = (TotalMoisture(),)
prognostic_vars(::NonEquilMoist) =
    (TotalMoisture(), LiquidMoisture(), IceMoisture())

prognostic_vars(::NoPrecipitation) = ()
prognostic_vars(::RainModel) = (Rain(),)
prognostic_vars(::RainSnowModel) = (Rain(), Snow())

prognostic_vars(::NoTracers) = ()
prognostic_vars(::NTracers{N}) where {N} = (Tracers{N}(),)

prognostic_vars(m::AtmosModel) = (
    Mass(),
    Momentum(),
    Energy(),
    prognostic_vars(m.moisture)...,
    prognostic_vars(m.precipitation)...,
    prognostic_vars(m.tracers)...,
    prognostic_vars(m.turbconv)...,
)

get_prog_state(state, ::Mass) = (state, :ρ)
get_prog_state(state, ::Momentum) = (state, :ρu)
get_prog_state(state, ::Energy) = (state, :ρe)
get_prog_state(state, ::TotalMoisture) = (state.moisture, :ρq_tot)
get_prog_state(state, ::Rain) = (state.precipitation, :ρq_rai)
get_prog_state(state, ::Snow) = (state.precipitation, :ρq_sno)
get_prog_state(state, ::Tracers{N}) where {N} = (state.tracers, :ρχ)
