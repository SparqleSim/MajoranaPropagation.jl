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
    sites_pattern = [searchsortedfirst(distinct_sites, s) for s in site_inds]
    raw_columns, n_fermions_canonical = _build_raw_columns(symbols, sites_acted_on, sites_pattern, is_spinful, theta)
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

# ============================================================================ #
#  Angle-free (reusable) lookup: build the transfer table once, evaluate at many
#  angles. The angle dependence of each coefficient is stored as a surrogate node
#  graph (`SurrogateCoeff`); `evaluate(table, theta)` materializes the concrete
#  numeric `FermionicRotationLookup` above, which is then propagated as usual.
# ============================================================================ #

"""
    SymbolicFermionicRotationLookup

Angle-free counterpart of [`FermionicRotationLookup`](@ref). Its transfer maps store, per output
term, a [`SurrogateCoeff`](@ref) carrying the symbolic angle dependence instead of a baked number.
Built once with the no-`theta` `FermionicRotationLookup(...)` constructors, then turned into a
concrete `FermionicRotationLookup` for a given angle with [`evaluate`](@ref).

`nparams` is the number of angles to supply at evaluation time (one per `symbol`; equal to 1 for
the single-symbol constructors).
"""
struct SymbolicFermionicRotationLookup{TM<:MajoranaTransferMap,STM<:MajoranaTransferMap}
    transfer_map::TM
    shifted_transfer_map::STM
    site_inds::Vector{Int}
    gate_modes::Vector{Int}
    is_spinful::Bool
    sites_acted_on::Int
    nparams::Int
end

function Base.show(io::IO, sym::SymbolicFermionicRotationLookup)
    print(io, "SymbolicFermionicRotationLookup(sites=$(sym.site_inds), is_spinful=$(sym.is_spinful), $(sym.transfer_map))")
end

"""
    FermionicRotationLookup(symbol, sites_acted_on, is_spinful, site_inds)
    FermionicRotationLookup(gate::FermionicRotation, is_spinful)

Angle-free constructors: build a reusable [`SymbolicFermionicRotationLookup`](@ref) (no `theta`).
Identical in structure to the `theta`-baking constructors above, but each basis string is
propagated through the surrogate so the table can later be [`evaluate`](@ref)d at many angles
without rebuilding it.
"""
function FermionicRotationLookup(symbols::Vector{Symbol}, sites_acted_on::Integer, is_spinful::Bool,
    site_inds::Vector{Int})

    distinct_sites = sort(unique(site_inds))
    gate_modes = _gate_modes(distinct_sites, is_spinful)
    sites_pattern = [searchsortedfirst(distinct_sites, s) for s in site_inds]

    raw_columns, n_fermions_canonical, nparams = _build_symbolic_columns(symbols, sites_acted_on, sites_pattern, is_spinful)
    TT_canonical = getinttype(n_fermions_canonical)

    transfer_map = MajoranaTransferMap(
        _finalize_symbolic_columns(raw_columns, TT_canonical, c -> TT_canonical(c), n_fermions_canonical, TT_canonical)
    )

    max_site = maximum(site_inds)
    n_fermions_shifted = is_spinful ? 2 * max_site : max_site
    TT_shifted = getinttype(n_fermions_shifted)
    shifted_transfer_map = MajoranaTransferMap(
        _finalize_symbolic_columns(raw_columns, TT_canonical, c -> _expand(c, gate_modes, TT_shifted), n_fermions_shifted, TT_shifted)
    )

    return SymbolicFermionicRotationLookup(transfer_map, shifted_transfer_map, site_inds, gate_modes, is_spinful, sites_acted_on, nparams)
end

function FermionicRotationLookup(symbol::Symbol, sites_acted_on::Integer, is_spinful::Bool, site_inds)
    return FermionicRotationLookup([symbol], sites_acted_on, is_spinful, collect(Int, site_inds))
end

function FermionicRotationLookup(gate::FermionicRotation, is_spinful::Bool)
    return FermionicRotationLookup(gate.symbol, length(gate.sites), is_spinful, gate.sites)
end

# Materialize a symbolic transfer map into a numeric one, reading the (already evaluated) node
# values: ctilde = a * inv_mu, dropping per-angle-zero entries to match the baked numeric table.
function _materialize_transfer_map(sym_map::MajoranaTransferMap{TT,SurrogateCoeff}) where {TT}
    ncol = ncolumns(sym_map)
    columns = Vector{Vector{Tuple{TT,ComplexF64}}}(undef, ncol)
    for col in 0:(ncol - 1)
        column = Tuple{TT,ComplexF64}[]
        for (cumulative, sc) in sym_map[col]
            ctilde = sc.path.node.cummulative_value * sc.inv_mu
            if abs(ctilde) < 1e-12
                continue
            end
            push!(column, (cumulative, ctilde))
        end
        columns[col + 1] = column
    end
    return MajoranaTransferMap(columns)
end

"""
    evaluate(sym::SymbolicFermionicRotationLookup, theta::Real)

Materialize the reusable symbolic lookup at rotation angle `theta`, returning a concrete numeric
[`FermionicRotationLookup`](@ref) (a `StaticGate`) ready to `propagate`. The node graphs are
evaluated once (shared between the canonical and shifted maps); the result is identical to
`FermionicRotationLookup(symbol, sites_acted_on, is_spinful, site_inds, theta)`.
"""
function evaluate(sym::SymbolicFermionicRotationLookup, theta::Real)
    thetas = fill(Float64(theta), sym.nparams)

    # evaluate the (shared) node graphs once; both maps reference the same node objects
    evaluate!([sc.path for (_, sc) in sym.transfer_map.entries], thetas)

    transfer_map = _materialize_transfer_map(sym.transfer_map)
    shifted_transfer_map = _materialize_transfer_map(sym.shifted_transfer_map)

    return FermionicRotationLookup(transfer_map, shifted_transfer_map, sym.site_inds, sym.gate_modes, sym.is_spinful, sym.sites_acted_on)
end

# ============================================================================ #
#  Canonical (placement-free) reusable lookup.
#
#  Unlike `SymbolicFermionicRotationLookup`, this stores ONLY the canonical node
#  columns (built on sites 1:k, independent of where the gate is placed). One table
#  serves every placement: `evaluate(table, theta, site_inds)` shifts the canonical
#  columns onto the actual modes and materializes a concrete `FermionicRotationLookup`.
# ============================================================================ #

"""
    CanonicalFermionicRotationLookup

Placement-free counterpart of [`SymbolicFermionicRotationLookup`](@ref): the symbolic transfer
table is built once on the canonical system, keyed by the *repetition pattern* of the gate's site
arguments (its `canonical_pattern`, e.g. `[1,1,2,3]` for a 4-argument gate whose first two arguments
coincide, or `[1,2]` for a plain 2-site gate). It does not commit to which actual sites the gate
acts on. Build it with the no-`theta` `FermionicRotationLookup(symbol, sites_acted_on, is_spinful;
sites_pattern=...)` constructor, then place + materialize it at an angle with
[`evaluate(::CanonicalFermionicRotationLookup, theta, site_inds)`](@ref).

`site_inds` supplied at evaluation must have a repetition pattern equal to the table's, i.e.
`canonical_pattern(site_inds) == table.canonical_pattern`; cache one table per pattern (e.g. keyed by
`canonical_pattern(site_inds)`). This supports repeated indices (`:four_body_term`, on-site terms)
and any ordering, including direction-sensitive gates.

Fields: `canonical_pattern` (gate args in canonical distinct-site labels) and `sites_acted_on`
(the number of *distinct* sites, `= maximum(canonical_pattern)`).
"""
struct CanonicalFermionicRotationLookup{TT<:Integer}
    raw_columns::Vector{Vector{Tuple{TT,MajoranaNodePathProperties}}}
    n_fermions_canonical::Int
    nparams::Int
    canonical_pattern::Vector{Int}
    sites_acted_on::Int
    is_spinful::Bool
end

function Base.show(io::IO, t::CanonicalFermionicRotationLookup)
    print(io, "CanonicalFermionicRotationLookup(pattern=$(t.canonical_pattern), is_spinful=$(t.is_spinful), $(length(t.raw_columns)) columns)")
end

"""
    canonical_pattern(site_inds)

Repetition pattern of `site_inds` (each site replaced by its rank among the distinct sites,
smallest = 1): `[5,5,7,8] -> [1,1,2,3]`, `[8,3] -> [2,1]`. Use it as the cache key for
[`CanonicalFermionicRotationLookup`](@ref) tables and to check [`evaluate`](@ref) compatibility.
"""
canonical_pattern(site_inds) = _rank_pattern(site_inds)

"""
    CanonicalFermionicRotationLookup(symbol, sites_acted_on, is_spinful; sites_pattern=nothing)

Build the canonical (placement-free) symbolic table. By default the gate acts on `sites_acted_on`
*distinct* sites (`canonical_pattern = 1:sites_acted_on`). For gates with repeated indices, pass a
representative `sites_pattern` (any indices — e.g. the actual `site_inds` like `[5,5,7,8]`, or a
bare pattern like `[1,1,2,3]`); only its repetition/order pattern is used (the literal values are
irrelevant), and `sites_acted_on` must equal its number of distinct sites. Also available via
`FermionicRotationLookup(symbol, sites_acted_on, is_spinful; sites_pattern=...)`.
"""
function CanonicalFermionicRotationLookup(symbols::Vector{Symbol}, sites_acted_on::Integer, is_spinful::Bool; sites_pattern=nothing)
    if isnothing(sites_pattern)
        pattern = collect(1:Int(sites_acted_on))
    else
        pattern = _rank_pattern(collect(Int, sites_pattern))
        maximum(pattern) == sites_acted_on || throw(ArgumentError(
            "sites_acted_on=$(sites_acted_on) must equal the number of distinct sites in " *
            "sites_pattern=$(sites_pattern) (= $(maximum(pattern)))."))
    end
    raw_columns, n_fermions_canonical, nparams = _build_symbolic_columns(symbols, sites_acted_on, pattern, is_spinful)
    return CanonicalFermionicRotationLookup(raw_columns, n_fermions_canonical, nparams, pattern, Int(sites_acted_on), is_spinful)
end

function CanonicalFermionicRotationLookup(symbol::Symbol, sites_acted_on::Integer, is_spinful::Bool; sites_pattern=nothing)
    return CanonicalFermionicRotationLookup([symbol], sites_acted_on, is_spinful; sites_pattern=sites_pattern)
end

# `FermionicRotationLookup(...)` with no `site_inds` and no `theta` => placement-free canonical table.
function FermionicRotationLookup(symbols::Vector{Symbol}, sites_acted_on::Integer, is_spinful::Bool; sites_pattern=nothing)
    return CanonicalFermionicRotationLookup(symbols, sites_acted_on, is_spinful; sites_pattern=sites_pattern)
end

function FermionicRotationLookup(symbol::Symbol, sites_acted_on::Integer, is_spinful::Bool; sites_pattern=nothing)
    return CanonicalFermionicRotationLookup(symbol, sites_acted_on, is_spinful; sites_pattern=sites_pattern)
end

"""
    evaluate(table::CanonicalFermionicRotationLookup, theta::Real, site_inds)

Place the canonical table on the actual `site_inds` and materialize a concrete numeric
[`FermionicRotationLookup`](@ref) at angle `theta`. `site_inds` must have the same repetition
pattern as the table (`canonical_pattern(site_inds) == table.canonical_pattern`); the result is then
identical to `FermionicRotation(symbol, site_inds)`.

The pattern (which arguments coincide, and the argument order) is baked into the table, while the
placement onto the actual modes is monotonic (canonical distinct site `j` → the `j`-th smallest
actual distinct site), which keeps the cumulative-string signs valid. This supports repeated indices
(e.g. `:four_body_term`, on-site terms) and any ordering, including direction-sensitive gates —
build one table per pattern, keyed by `canonical_pattern(site_inds)`.
"""
function evaluate(table::CanonicalFermionicRotationLookup{TT}, theta::Real, site_inds) where {TT}
    sites = collect(Int, site_inds)
    n_args = length(table.canonical_pattern)
    if length(sites) != n_args
        throw(ArgumentError("evaluate expects $n_args site index/indices to match the table's pattern $(table.canonical_pattern), got $(length(sites)): $sites"))
    end
    if _rank_pattern(sites) != table.canonical_pattern
        throw(ArgumentError(
            "site_inds $sites has repetition pattern $(_rank_pattern(sites)), which does not match the " *
            "table's canonical_pattern $(table.canonical_pattern). Build/select a table whose pattern " *
            "equals `canonical_pattern(site_inds)`."))
    end

    thetas = fill(Float64(theta), table.nparams)
    # monotonic placement: canonical distinct site j -> j-th smallest actual distinct site
    distinct_sites = sort(unique(sites))
    gate_modes = _gate_modes(distinct_sites, table.is_spinful)

    n_fermions_shifted = table.is_spinful ? 2 * maximum(sites) : maximum(sites)
    TT_shifted = getinttype(n_fermions_shifted)

    canonical_sym = MajoranaTransferMap(
        _finalize_symbolic_columns(table.raw_columns, TT, c -> TT(c), table.n_fermions_canonical, TT)
    )
    shifted_sym = MajoranaTransferMap(
        _finalize_symbolic_columns(table.raw_columns, TT, c -> _expand(c, gate_modes, TT_shifted), n_fermions_shifted, TT_shifted)
    )

    # evaluate the (shared) node graphs once; both maps reference the same node objects
    evaluate!([sc.path for (_, sc) in canonical_sym.entries], thetas)

    transfer_map = _materialize_transfer_map(canonical_sym)
    shifted_transfer_map = _materialize_transfer_map(shifted_sym)

    return FermionicRotationLookup(transfer_map, shifted_transfer_map, sites, gate_modes, table.is_spinful, n_args)
end
