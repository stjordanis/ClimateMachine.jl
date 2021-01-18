##### Sum wrapper

using StaticArrays

export Σfluxes, Σsources

"""
    flux

An individual flux.
See [`BalanceLaw`](@ref) for more info.
"""
function flux end

"""
    source

An individual source.
See [`BalanceLaw`](@ref) for more info.
"""
function source end

"""
    ntuple_sum(nt::NTuple{N,T}) where {N, T}

sum of `NTuple`, which requires more strict
type input than `sum`. This is added to better
synchronize the success/failure between CPU/GPU
runs to help improve debugging.
"""
ntuple_sum(nt::NTuple{N, T}) where {N, T} = sum(nt)

"""
    Σfluxes(fluxes::NTuple, args...)

Sum of the fluxes where
 - `fluxes` is an `NTuple{N, TendencyDef{Flux{O}, PV}} where {N, PV, O}`
"""
function Σfluxes(
    pv::PV,
    fluxes::NTuple{N, TendencyDef{Flux{O}, PV}},
    args...,
) where {N, PV, O}
    return ntuple_sum(ntuple(Val(N)) do i
        flux(fluxes[i], args...)
    end)
end

# TODO: is there a cleaner way?
Σfluxes(
    pv::PV,
    fluxes::NTuple{0, TendencyDef{Flux{O}, PV}},
    args...,
) where {PV, O} = Σfluxes(Val(n_components(pv)), fluxes, args...)

Σfluxes(
    ::Val{1},
    fluxes::NTuple{0, TendencyDef{Flux{O}, PV}},
    args...,
) where {PV, O} = SArray{Tuple{3}}(ntuple(i -> 0, 3))

Σfluxes(
    ::Val{nc},
    fluxes::NTuple{0, TendencyDef{Flux{O}, PV}},
    args...,
) where {nc, PV, O} = SArray{Tuple{nc}}(ntuple(i -> 0, nc))

"""
    Σsources(sources::NTuple, args...)

Sum of the sources where
 - `sources` is an `NTuple{N, TendencyDef{Source, PV}} where {N, PV}`
"""
function Σsources(
    pv::PV,
    sources::NTuple{N, TendencyDef{Source, PV}},
    args...,
) where {N, PV}
    return ntuple_sum(ntuple(Val(N)) do i
        source(sources[i], args...)
    end)
end
Σsources(
    pv::PV,
    sources::NTuple{0, TendencyDef{Source, PV}},
    args...,
) where {PV} = Σsources(Val(n_components(pv)), sources, args...)

Σsources(
    ::Val{1},
    sources::NTuple{0, TendencyDef{Source, PV}},
    args...,
) where {PV} = 0

Σsources(
    ::Val{nc},
    sources::NTuple{0, TendencyDef{Source, PV}},
    args...,
) where {nc, PV} = SArray{Tuple{nc}}(ntuple(i -> 0, nc))
