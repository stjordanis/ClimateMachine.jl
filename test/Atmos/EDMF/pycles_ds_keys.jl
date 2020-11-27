var_map(s::String) = var_map(Val(Symbol(s)))
var_map(::Val{T}) where {T} = nothing

var_map(::Val{Symbol("ρ")})                                = ("rho", ())
var_map(::Val{Symbol("ρu[1]")})                            = ("u_mean", (:ρ,))
var_map(::Val{Symbol("ρu[2]")})                            = ("v_mean", (:ρ,))
# var_map(::Val{Symbol("aux.moisture.θ_liq")})               = ("thetali_mean", (:ρ,:a))
var_map(::Val{Symbol("moisture.ρq_tot")})                  = ("qt_mean", (:ρ,))
var_map(::Val{Symbol("turbconv.updraft[1].ρa")})           = ("updraft_fraction", (:ρ,))
var_map(::Val{Symbol("turbconv.updraft[1].ρaw")})          = ("updraft_w", (:ρ,:a))
var_map(::Val{Symbol("turbconv.updraft[1].ρaq_tot")})      = ("updraft_qt", (:ρ,:a))
var_map(::Val{Symbol("turbconv.updraft[1].ρaθ_liq")})      = ("updraft_thetali", (:ρ,:a))
var_map(::Val{Symbol("turbconv.environment.ρatke")})       = ("tke_mean", (:ρ,:a))
var_map(::Val{Symbol("turbconv.environment.ρaθ_liq_cv")})  = ("env_thetali2", (:ρ,:a))
var_map(::Val{Symbol("turbconv.environment.ρaq_tot_cv")})  = ("env_qt2", (:ρ,:a))

# PyCLES dataset keys:
# "rho"
# "p0"
# "Hvar_mean"
# "QTvar_mean"
# "env_Hvar"
# "env_QTvar"
# "env_HQTcov"
# "W_third_m"
# "H_third_m"
# "QT_third_m"
# "massflux"
# "massflux_h"
# "massflux_qt"
# "total_flux_h"
# "total_flux_qt"
# "diffusive_flux_h"
# "diffusive_flux_qt"
# "total_flux_u"
# "total_flux_v"
# "diffusive_flux_u"
# "diffusive_flux_v"
# "massflux_u"
# "massflux_v"
# "buoyancy_mean"
# "cloud_fraction"
# "temperature_mean"
# "updraft_ddz_p_alpha"
# "thetali_mean"
# "qt_mean"
# "ql_mean"
# "u_mean"
# "v_mean"
# "tke_mean"
# "v_translational_mean"
# "u_translational_mean"
# "updraft_buoyancy"
# "updraft_fraction"
# "env_thetali"
# "updraft_thetali"
# "env_qt"
# "updraft_qt"
# "env_ql"
# "env_buoyancy"
# "updraft_ql"
# "qr_mean"
# "env_qr"
# "updraft_qr"
# "updraft_w"
# "updraft_RH"
# "env_RH"
# "env_w"
# "thetali_mean2"
# "qt_mean2"
# "env_thetali2"
# "env_qt2"
# "env_qt_thetali"
# "tke_prod_A"
# "tke_prod_B"
# "tke_prod_D"
# "tke_prod_P"
# "tke_prod_T"
# "tke_prod_S"
# "tke_nd_mean"