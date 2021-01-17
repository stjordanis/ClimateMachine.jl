# include("src\\Diagnostics\\tendency_diagnostics.jl")
# include("test\\Atmos\\EDMF\\bomex_edmf.jl")
# let

tendency_diagnostics(solver_config, diag_in) = tendency_diagnostics(
    solver_config.dg.balance_law,
    diag_in.prog,
    diag_in.aux,
    diag_in.∇flux,
    diag_in.hyperdiff,
    400.0,
    VerticalDirection(),
)

# Process prognostic variable name:
params(::T) where {T} = T.parameters
unval(::Val{i}) where {i} = i
param_suffix(p, ::Val{n}) where {n} = Symbol(:_, string(unval(p[1])))
param_suffix(p, ::Val{0}) = Symbol()
param_suffix(p) = param_suffix(p, Val(length(p)))
prog_name(pv::PV) where {PV} = Symbol(nameof(PV), param_suffix(params(pv)))
tend_name(t) = nameof(typeof(t))

"""
    tend_vals(bl, tend_fun, tend_type, args)

A nested NamedTuple of tendency terms
per prognostic variable.
"""
tend_vals(bl, tend_fun, tend_type, args) = [
    begin
        eqts = eq_tends(pv, bl, tend_type)
        source_term_keys = tend_name.(eqts)
        source_term_vals = ntuple(Val(length(eqts))) do i
            tend_fun(eqts[i], bl, args)
        end
        (; zip(source_term_keys, source_term_vals)...)
    end for pv in prognostic_vars(bl)
]

"""
    tendency_diagnostics

A flattened NamedTuple of tendency terms.
The keys use the following convention:
 - `key1_key2_key3[_key4]`
where
 - `key1`: tendency type
 - `key2`: prognostic variable
 - `key3`: tendency term
 - `key4`: vector, or flux, component (if applicable)

# Example:
```julia
diag = tendency_diagnostics(...)

# Gravity source for x-component momentum:
diag.src_Momentum_Gravity_1
```

TODO: Add options for filtering / selecting
      specific terms

!!! warn
    This function returns a _very_ large NamedTuple,
    is *extremely costly*, and *should be used
    infrequently while debugging small problems*.
    For example, when first adding this functionality,
    in `bomex_edmf.jl`, this NamedTuple had 254 keys
    and values.
"""
function tendency_diagnostics end

tendency_diagnostics(solver_config; kwargs...) = tendency_diagnostics(
    solver_config.dg.grid,
    solver_config.dg.balance_law,
    gettime(solver_config.solver),
    solver_config.dg.direction;
    prognostic = solver_config.Q,
    auxiliary = solver_config.dg.state_auxiliary,
    diffusive = solver_config.dg.state_gradient_flux,
    hyperdiffusive = solver_config.dg.states_higher_order[2],
    kwargs...,
)

function tendency_diagnostics(
    grid::DiscontinuousSpectralElementGrid,
    bl::BalanceLaw,
    t::Real,
    direction;
    kwargs...,
)
    return [
        begin
            tendency_diagnostics(
                bl,
                state,
                aux,
                diffusive,
                hyperdiffusive,
                t,
                direction,
            )
        end for local_states in NodalStack(bl, grid; kwargs...)
    ]
end

function tendency_diagnostics(
    bl::BalanceLaw,
    state,
    aux,
    diffusive,
    hyperdiffusive,
    t::Real,
    direction,
)
    _args_fx1 = (; state, aux, t, direction)
    _args_fx2 = (; state, aux, t, diffusive, hyperdiffusive)
    _args_src = (; state, aux, t, direction, diffusive)

    cache_fx1 = precompute(bl, _args_fx1, Flux{FirstOrder}())
    cache_fx2 = precompute(bl, _args_fx2, Flux{SecondOrder}())
    cache_src = precompute(bl, _args_src, Source())

    # cache_fx1, cache_fx2, and cache_src have overlapping
    # data, only need one copy, so merge:
    cache = merge(cache_fx1, cache_fx2, cache_src)
    z = altitude(bl, aux)

    # ----------- Compute tendencies
    args_fx1 = merge(_args_fx1, (; precomputed = cache_fx1))
    args_fx2 = merge(_args_fx2, (; precomputed = cache_fx2))
    args_src = merge(_args_src, (; precomputed = cache_src))

    prog_vars = prog_name.(prognostic_vars(bl))

    flux_O1_vals = tend_vals(bl, flux, Flux{FirstOrder}(), args_fx1)
    flux_O2_vals = tend_vals(bl, flux, Flux{SecondOrder}(), args_fx2)
    source_vals = tend_vals(bl, source, Source(), args_src)

    flx1 = (; zip(prog_vars, flux_O1_vals)...)
    flx2 = (; zip(prog_vars, flux_O2_vals)...)
    src = (; zip(prog_vars, source_vals)...)

    nt = (;
        z = altitude(bl, aux),
        prog = flattened_named_tuple(state), # Vars -> flattened NamedTuples
        aux = flattened_named_tuple(aux), # Vars -> flattened NamedTuples
        ∇flux = flattened_named_tuple(diffusive), # Vars -> flattened NamedTuples
        hyperdiff = flattened_named_tuple(hyperdiffusive), # Vars -> flattened NamedTuples
        cache = cache,
        flx1 = flx1,
        flx2 = flx2,
        src = src,
    )

    # Flatten top level:
    return flattened_named_tuple(nt)
end


# diag_vs_z = single_stack_diagnostics(solver_config)[1]

# using Test
# @test diag_vs_z.hyperdiff == nothing
# @test diag_vs_z.cache_ts isa ThermodynamicState
# @test diag_vs_z.cache_turbulence_τ_9 isa AbstractFloat

# diag_in = diagnostics_input(solver_config)
# diag_vs_z = tendency_diagnostics(solver_config, diag_in)

# for k in keys(diag_vs_z)
#     @show k
# end

# nothing
# end
