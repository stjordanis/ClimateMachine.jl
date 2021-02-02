#### Surface model kernels

using Statistics
using ClimateMachine.SurfaceFluxes
const ArrayType = ClimateMachine.array_type()

"""
    subdomain_surface_values(
        surf::SurfaceModel,
        turbconv::EDMF{FT},
        atmos::AtmosModel{FT},
        state::Vars,
        aux::Vars,
        state_int::Vars,
        aux_int::Vars,
        zLL::FT,
    ) where {FT}

Returns the surface values of updraft area fraction, updraft
liquid water potential temperature (`θ_liq`), updraft total
water specific humidity (`q_tot`), environmental variances of
`θ_liq` and `q_tot`, environmental covariance of `θ_liq` with
`q_tot`, and environmental TKE, given:
 - `surf`, a `SurfaceModel`
 - `turbconv`, an `EDMF` model
 - `atmos`, an `AtmosModel`
 - `state`, state variables
 - `aux`, auxiliary variables
 - `state_int`, state variables at first interior point
 - `aux_int`, auxiliary variables at first interior point
 - `zLL`, height of the lowest nodal level
"""
function subdomain_surface_values(
    surf::SurfaceModel,
    turbconv::EDMF{FT},
    atmos::AtmosModel{FT},
    state::Vars,
    aux::Vars,
    state_int::Vars,
    aux_int::Vars,
    zLL::FT,
) where {FT}

    turbconv = atmos.turbconv
    N_up = n_updrafts(turbconv)
    gm = state
    gm_int = state_int
    # TODO: change to new_thermo_state
    ts = recover_thermo_state(atmos, state, aux)
    ts_int = recover_thermo_state(atmos, state_int, aux_int)
    q = PhasePartition(ts)
    _cp_m = cp_m(ts)
    lv = latent_heat_vapor(ts)
    Π = exner(ts)
    ρ_inv = 1 / gm.ρ
    ρ_int_inv = 1 / gm.ρ

    surface_scalar_coeff = turbconv.surface.scalar_coeff

    # Retrieve surface fluxes based on MO similarity
    ## Initial guesses for MO parameters
    LMO_init = eps(FT)
    u_star_init = FT(0.1)
    th_star_init = eps(FT)
    qt_star_init = eps(FT)
    x_init = ArrayType(FT[LMO_init, u_star_init, th_star_init, qt_star_init])

    # Surface values for variables: To be revised
    u_sfc = FT(0)
    thv_sfc = virtual_pottemp(ts)
    qt_sfc = total_specific_humidity(ts)
    x_s = ArrayType(FT[u_sfc, thv_sfc, qt_sfc])
    VDSE_scale = _cp_m * virtual_pottemp(ts) # Revise

    # Avg values in first cell for variables: To be revised
    qt_ave = total_specific_humidity(ts_int)
    vdse_ave = _cp_m * virtual_pottemp(ts_int)
    u_ave = sqrt(
        (gm_int.ρu[1] * ρ_int_inv) * (gm_int.ρu[1] * ρ_int_inv) +
        (gm_int.ρu[2] * ρ_int_inv) * (gm_int.ρu[2] * ρ_int_inv),
    )
    x_ave = ArrayType(FT[u_ave, vdse_ave, qt_ave])

    z_rough = ArrayType(FT[surf.z_0, surf.z_0, surf.z_0])

    surf_flux_cond = surface_conditions(
        atmos.param_set,
        x_init,
        x_ave,
        x_s,
        z_rough,
        VDSE_scale,
        zLL,
        DGScheme(),
    )

    θ_liq_surface_flux = surf_flux_cond.flux[2]
    q_tot_surface_flux = surf_flux_cond.flux[3]

    oblength = surf_flux_cond.L_MO
    ustar = surf_flux_cond.x_star[1]
    println("Obukhov length is", oblength)
    println("Friction velocity is", ustar)

    unstable = oblength < 0
    fact = unstable ? (1 - surf.ψϕ_stab * zLL / oblength)^(-FT(2 // 3)) : 1
    tke_fact = unstable ? cbrt(zLL / oblength * zLL / oblength) : 0
    ustar² = ustar^2
    θ_liq_cv = 4 * (θ_liq_surface_flux * θ_liq_surface_flux) / (ustar²) * fact
    q_tot_cv = 4 * (q_tot_surface_flux * q_tot_surface_flux) / (ustar²) * fact
    θ_liq_q_tot_cv =
        4 * (θ_liq_surface_flux * q_tot_surface_flux) / (ustar²) * fact
    tke = ustar² * (surf.κ_star² + tke_fact)

    a_up_surf = ntuple(i -> FT(surf.a / N_up), N_up)
    e_int = internal_energy(atmos, state, aux)
    ts_new = new_thermo_state(atmos, state, aux)
    θ_liq = liquid_ice_pottemp(ts_new)

    upd_θ_liq_surf = ntuple(N_up) do i
        θ_liq + surface_scalar_coeff[i] * sqrt(max(θ_liq_cv, 0))
    end

    ρq_tot = atmos.moisture isa DryModel ? FT(0) : gm.moisture.ρq_tot
    q_tot_up_surf = ntuple(N_up) do i
        ρq_tot * ρ_inv + surface_scalar_coeff[i] * sqrt(max(q_tot_cv, 0))
    end

    return (
        a_up_surf = a_up_surf,
        upd_θ_liq_surf = upd_θ_liq_surf,
        q_tot_up_surf = q_tot_up_surf,
        θ_liq_cv = θ_liq_cv,
        q_tot_cv = q_tot_cv,
        θ_liq_q_tot_cv = θ_liq_q_tot_cv,
        tke = tke,
        ob_length = oblength,
        u_star = ustar,
    )
end;

"""
    percentile_bounds_mean_norm(
        low_percentile::FT,
        high_percentile::FT,
        n_samples::Int,
    ) where {FT <: AbstractFloat}

Returns the mean of all instances of a standard Gaussian random
variable that have a CDF higher than low_percentile and lower
than high_percentile, given a total of n_samples of the standard
Gaussian, given:
 - `low_percentile`, lower limit of the CDF
 - `high_percentile`, higher limit of the CDF
 - `n_samples`, the total number of samples drawn from the Gaussian
"""
function percentile_bounds_mean_norm(
    low_percentile::FT,
    high_percentile::FT,
    n_samples::Int,
) where {FT <: AbstractFloat}
    x = rand(Normal(), n_samples)
    xp_low = quantile(Normal(), low_percentile)
    xp_high = quantile(Normal(), high_percentile)
    filter!(y -> xp_low < y < xp_high, x)
    return Statistics.mean(x)
end
