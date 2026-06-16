"""
    FermionicRotationLookup(symbol::Symbol, sites_acted_on, is_spinful, site_inds, theta)

A `StaticGate` that applies a `FermionicRotation` via a precomputed transfer map (lookup
table), analogous to `PauliPropagation.TransferMapGate`. The rotation angle `theta` is baked
into the table at construction time.

The gate acts on `sites_acted_on` lattice sites placed at `site_inds` of the full system
(spinful or spinless according to `is_spinful`). It is intended as a drop-in alternative to
applying `FermionicRotation(symbol, site_inds)` with parameter `theta`, and produces the same
Majorana sum.

Fields
- `transfer_map`: the canonical transfer map (cumulative strings on canonical modes `1:k`).
- `shifted_transfer_map`: the transfer map with cumulative strings shifted onto `site_inds`.
- `site_inds`: the sites the gate acts on, in the order given (order matters for e.g. :hopupdn).
- `gate_modes`: the (sorted, 0-based) Majorana mode bit positions of the gate.
- `is_spinful`, `sites_acted_on`: gate metadata.
"""
struct FermionicRotationLookup{TM<:MajoranaTransferMap,STM<:MajoranaTransferMap} <: StaticGate
    transfer_map::TM
    shifted_transfer_map::STM
    site_inds::Vector{Int}
    gate_modes::Vector{Int}
    is_spinful::Bool
    sites_acted_on::Int
end

function FermionicRotationLookup(symbols::Vector{Symbol}, sites_acted_on::Integer, is_spinful::Bool,
    site_inds::Vector{Int}, theta::Real)

    # `site_inds` may any collection of sites so that `FermionicRotation(symbol, site_inds)`
    # is a valid operator with that repetition)
    distinct_sites = sort(unique(site_inds))

    gate_modes = _gate_modes(distinct_sites, is_spinful)

    # Build the lookup table on a compact canonical system: each actual site is mapped to its rank
    # within the distinct sites (1 = smallest, ..., k = largest). This compacts the system to sites
    # 1:k while preserving both the order and any repetition of `site_inds`, e.g. [5, 2] -> [2, 1]
    # and [2, 2] -> [1, 1]. The order-preserving (monotonic) mapping keeps the canonical
    # coefficients valid once shifted onto the actual modes.
    canonical_sites = [searchsortedfirst(distinct_sites, s) for s in site_inds]
    raw_columns, n_fermions_canonical = _build_raw_columns(symbols, sites_acted_on, canonical_sites, is_spinful, theta)
    TT_canonical = getinttype(n_fermions_canonical)

    # canonical transfer map (cumulative strings on modes 1:k)
    transfer_map = MajoranaTransferMap(
        _finalize_columns(raw_columns, TT_canonical, c -> TT_canonical(c), n_fermions_canonical, TT_canonical)
    )

    # shifted transfer map (cumulative strings on the actual modes `site_inds`)
    max_site = maximum(site_inds)
    n_fermions_shifted = is_spinful ? 2 * max_site : max_site
    TT_shifted = getinttype(n_fermions_shifted)
    shifted_transfer_map = MajoranaTransferMap(
        _finalize_columns(raw_columns, TT_canonical, c -> _expand(c, gate_modes, TT_shifted), n_fermions_shifted, TT_shifted)
    )

    return FermionicRotationLookup(transfer_map, shifted_transfer_map, site_inds, gate_modes, is_spinful, sites_acted_on)
end

function FermionicRotationLookup(symbol::Symbol, sites_acted_on::Integer, is_spinful::Bool, site_inds, theta::Real)
    return FermionicRotationLookup([symbol], sites_acted_on, is_spinful, collect(Int, site_inds), theta)
end

"""
    FermionicRotationLookup(gate::FermionicRotation, n_sites, is_spinful, theta)

Convenience constructor that builds a lookup gate equivalent to `gate` acting with `theta`.
"""
function FermionicRotationLookup(gate::FermionicRotation, is_spinful::Bool, theta::Real)
    return FermionicRotationLookup(gate.symbol, length(gate.sites), is_spinful, gate.sites, theta)
end

"""
The conjugate action `U' O U` of the fermionic gate is read off the transfer map.

For each observable string `O` we extract the restricted string on the gate's modes (the
column index), and for each `(cumulative_string, ctilde)` entry of that column compute
`(η, O_new) = ms_mult(cumulative_string, O)`. The full-string-dependent sign `η` combines
with the stored prefactor `ctilde` to give the new coefficient `coeff * real(η * ctilde)`.
Strings with no support on the gate's modes (column 0) are left unchanged.
"""
function PropagationBase.applytoall!(gate::FermionicRotationLookup, prop_cache::MajoranaPropagationCache; kwargs...)
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    n_fermions = nfermions(msum)
    smap = gate.shifted_transfer_map
    gate_modes = gate.gate_modes

    for (ms_int, coeff) in msum
        column_index = _compress(ms_int, gate_modes)

        # no support on the gate's modes => the gate commutes through, string is unchanged
        if column_index == 0
            add!(aux_msum, ms_int, coeff)
            continue
        end

        TT = typeof(ms_int)
        for (cumulative, ctilde) in smap[column_index]
            sign, new_ms = ms_mult(convert(TT, cumulative), ms_int, n_fermions)
            add!(aux_msum, new_ms, coeff * real(sign * ctilde))
        end
    end

    # everything was moved into the auxiliary sum; mirror the default applytoall! bookkeeping
    empty!(msum)
    PropagationBase.swapsums!(prop_cache)

    return
end
