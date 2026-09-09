###
##
# Fused `applymergetruncate!` overloads for `VectorMajoranaPropagationCache`.
# A rotation is applied in one counting and one writing pass. Products that fail a truncation
# reading only the Majorana string (`max_weight`, `max_unpaired`) are dropped as they are made and
# never written or sorted; every truncation is then applied inside the tail merge that follows, so
# no separate `truncate!` pass runs. The results are identical to the default path.
# Only used when `fused=true`; otherwise falls through (via `invoke`) to the default methods.
##
###

# truncation kwargs that keep everything
const _NOTRUNC = (min_abs_coeff=0.0, max_weight=Inf, max_unpaired=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing)

"""
    applymergetruncate!(gate::MajoranaRotation, prop_cache::VectorMajoranaPropagationCache, theta; fused::Bool=false, kwargs...)

Fused, task-partitioned overload of `applymergetruncate!` for `MajoranaRotation` -- see file header.
`force_truncate=true` tests every term against the truncations, including those the gate leaves alone.
"""
function PropagationBase.applymergetruncate!(gate::MajoranaRotation, prop_cache::MP.VectorMajoranaPropagationCache, theta;
    fused::Bool=false,
    min_abs_coeff::Real=1e-10, max_weight::Real=Inf, max_unpaired::Real=Inf, unpaired_mask=nothing,
    max_freq::Real=Inf, max_sins::Real=Inf, customtruncfunc=nothing,
    thread::Bool=true, force_truncate::Bool=false, truncate_in_merge::Bool=true, kwargs...)

    if !fused
        return invoke(PropagationBase.applymergetruncate!,
            Tuple{MajoranaRotation,MP.AbstractMajoranaPropagationCache,typeof(theta)},
            gate, prop_cache, theta;
            min_abs_coeff, max_weight, max_unpaired, unpaired_mask, max_freq, max_sins, customtruncfunc, thread, kwargs...)
    end

    _fusedrotation!(gate.ms_int, prop_cache, cos(theta), sin(theta), Val(false);
        min_abs_coeff, max_weight, max_unpaired, unpaired_mask, max_freq, max_sins, customtruncfunc,
        thread, force_truncate, truncate_in_merge)

    return prop_cache
end

"""
    applymergetruncate!(gate::ImaginaryMajoranaRotation, prop_cache::VectorMajoranaPropagationCache, beta; fused::Bool=false, kwargs...)

Fused, task-partitioned overload of `applymergetruncate!` for `ImaginaryMajoranaRotation`, sharing its kernel with the `MajoranaRotation` overload.
"""
function PropagationBase.applymergetruncate!(gate::ImaginaryMajoranaRotation, prop_cache::MP.VectorMajoranaPropagationCache, beta;
    fused::Bool=false,
    min_abs_coeff::Real=1e-10, max_weight::Real=Inf, max_unpaired::Real=Inf, unpaired_mask=nothing,
    max_freq::Real=Inf, max_sins::Real=Inf, customtruncfunc=nothing,
    thread::Bool=true, force_truncate::Bool=false, truncate_in_merge::Bool=true, kwargs...)

    if !fused
        # there is no Majorana-specific default method for a bare imaginary rotation
        return invoke(PropagationBase.applymergetruncate!,
            Tuple{ImaginaryMajoranaRotation,AbstractPropagationCache,typeof(beta)},
            gate, prop_cache, beta;
            min_abs_coeff, max_weight, max_unpaired, unpaired_mask, max_freq, max_sins, customtruncfunc, thread, kwargs...)
    end

    _fusedrotation!(gate.ms_int, prop_cache, cosh(beta), -sinh(beta), Val(true);
        min_abs_coeff, max_weight, max_unpaired, unpaired_mask, max_freq, max_sins, customtruncfunc,
        thread, force_truncate, truncate_in_merge)

    return prop_cache
end

### Shared kernel

# Applies the even-weight rotation string `gate_ms` to every active term and merges the products
# back in. `Commuting` selects the terms that branch: those anticommuting with the gate for
# real-time rotations (`false`), those commuting with it for imaginary-time ones (`true`). A
# branching term keeps its string and has its coefficient scaled by `kept_val` in place; its
# product gets `new_val * sign` and is appended, unless it fails a term-only truncation.
# `truncate_in_merge=false` merges without truncating, for callers that truncate later themselves
# (after renormalizing); the term-only truncations of the products do not depend on the scale.
# Returns `(n_branched, n_products)`.
function _fusedrotation!(gate_ms::TT, prop_cache::MP.VectorMajoranaPropagationCache, kept_val, new_val, ::Val{Commuting};
    min_abs_coeff::Real, max_weight::Real, max_unpaired::Real, unpaired_mask, max_freq::Real, max_sins::Real, customtruncfunc,
    thread::Bool, force_truncate::Bool=false, truncate_in_merge::Bool=true) where {TT,Commuting}

    n_old = activesize(prop_cache)
    if n_old == 0
        return (0, 0)
    end

    mask = MP._unpairedmask(prop_cache, unpaired_mask, max_unpaired)
    gate_ms_ps = MP.compute_parity_bits_and_shift(gate_ms, 2 * nfermions(prop_cache))
    omega_l_gate = MP.omega_L_mult(gate_ms)

    # the XOR-sorted tail merge needs the whole active range sorted before the apply
    sorted_before = sortedprefix(mainsum(prop_cache)) == n_old

    task_partitioner, n_tasks = PropagationBase._preparetasks(n_old, thread)
    branched_counts = Vector{Int}(undef, n_tasks)
    product_counts = Vector{Int}(undef, n_tasks)

    # dry run: each task counts the terms it branches and the products it will append
    _fusedpass!(prop_cache, task_partitioner, n_tasks, branched_counts, product_counts, zeros(Int, n_tasks + 1),
        gate_ms, gate_ms_ps, omega_l_gate, kept_val, new_val, max_weight, max_unpaired, mask, Val(Commuting), Val(false))

    n_branched = sum(branched_counts)
    if n_branched == 0 && !force_truncate
        return (0, 0)
    end

    new_offsets = PropagationBase._offsetsfromcounts(product_counts)
    n_products = new_offsets[end] - 1

    resize_factor = 1.5
    if capacity(prop_cache) < n_old + n_products
        resize!(prop_cache, round(Int, (n_old + n_products) * resize_factor))
    end

    # real pass: scale the branching coefficients where they lie and append each task's products
    # at its own offset past the end
    _fusedpass!(prop_cache, task_partitioner, n_tasks, branched_counts, product_counts, new_offsets,
        gate_ms, gate_ms_ps, omega_l_gate, kept_val, new_val, max_weight, max_unpaired, mask, Val(Commuting), Val(true))

    # the old terms kept their strings and their order, so the sorted prefix stands
    setactivesize!(prop_cache, n_old + n_products)

    truncfunc = truncate_in_merge ?
                MP._truncationfunc(mask; min_abs_coeff, max_weight, max_unpaired, max_freq, max_sins, customtruncfunc) :
                nothing
    MP._mergeafterapply!(prop_cache, sorted_before ? gate_ms : nothing; thread, truncfunc)

    # the merge returns right away when nothing was appended, skipping its test of the old terms
    if n_products == 0 && truncfunc !== nothing
        truncate!(truncfunc, prop_cache; thread)
    end

    return (n_branched, n_products)
end

# One task-partitioned pass over the active terms; only counts when `DoWrite` is false. Fetches
# the backing arrays itself, after any resize, so that the closure never captures a reassigned
# variable.
function _fusedpass!(prop_cache, task_partitioner, n_tasks, branched_counts, product_counts, new_offsets,
    gate_ms, gate_ms_ps, omega_l_gate, kept_val, new_val, max_weight, max_unpaired, mask, ::Val{Commuting}, ::Val{DoWrite}) where {Commuting,DoWrite}

    n_old = activesize(prop_cache)
    terms, coeffs, _, _ = PropagationBase._mainauxarrays(prop_cache)

    AK.itask_partition(n_tasks, n_tasks, 1) do task_id, _
        rng = task_partitioner[task_id]
        branched_counts[task_id], product_counts[task_id] = _fusedbranchwrite!(terms, coeffs, n_old + new_offsets[task_id], rng.start, rng.stop,
            gate_ms, gate_ms_ps, omega_l_gate, kept_val, new_val, max_weight, max_unpaired, mask, Val(Commuting), Val(DoWrite))
    end

    return
end

# Walks terms[lo:hi] and branches the terms `Commuting` selects. Returns `(n_branched, n_products)`.
@inline function _fusedbranchwrite!(terms, coeffs, new_start, lo, hi,
    gate_ms::TT, gate_ms_ps::TT, omega_l_gate::Int, kept_val, new_val, max_weight, max_unpaired, mask,
    ::Val{Commuting}, ::Val{DoWrite}) where {TT,Commuting,DoWrite}

    n_branched = 0
    new_pos = new_start

    @inbounds for ii in lo:hi
        term = terms[ii]
        MP._commutes_evengate(term, gate_ms) == Commuting || continue
        n_branched += 1

        coeff = coeffs[ii]
        DoWrite && (coeffs[ii] = _keptcoeff(coeff, kept_val, Val(Commuting)))

        new_term, sign = _gateproduct(term, gate_ms, gate_ms_ps, omega_l_gate, Val(Commuting))
        if !MP._truncateterm(new_term, max_weight, max_unpaired, mask)
            new_pos = PropagationBase._writeandadvance!(terms, coeffs, new_pos, new_term, _newcoeff(coeff, new_val * sign, Val(Commuting)), Val(DoWrite))
        end
    end

    return (n_branched, new_pos - new_start)
end

# real-time rotations go through `_applycos`/`_applysin`, which also advance the counters of a
# `MajoranaFrequencyTracker`; imaginary-time rotations scale plainly, as their default kernel does
@inline _gateproduct(term, gate_ms, gate_ms_ps, omega_l_gate, ::Val{false}) = MP._rotationproduct_evengate(term, gate_ms, gate_ms_ps, omega_l_gate)
@inline _gateproduct(term, gate_ms, gate_ms_ps, omega_l_gate, ::Val{true}) = MP._rotationproduct_evengate_commuting(term, gate_ms, gate_ms_ps, omega_l_gate)
@inline _keptcoeff(coeff, val, ::Val{false}) = MP._applycos(coeff, val)
@inline _keptcoeff(coeff, val, ::Val{true}) = coeff * val
@inline _newcoeff(coeff, val, ::Val{false}) = MP._applysin(coeff, val)
@inline _newcoeff(coeff, val, ::Val{true}) = coeff * val

### Fermionic gates

"""
    applymergetruncate!(gate::FermionicRotation, prop_cache, theta; fused::Bool=false, truncate_each_mr=nothing, kwargs...)

Fused overload for `FermionicRotation` on a `MajoranaSum` or `VectorMajoranaSum` cache: applies the Majorana rotations the gate decomposes into with the fused `MajoranaRotation` overload of the cache.
As in the default method, a gate whose rotations are only number-preserving together (hoppings) is truncated once after its last rotation, any other after every rotation.
"""
function PropagationBase.applymergetruncate!(gate::FermionicRotation, prop_cache::Union{MajoranaPropagationCache,MP.VectorMajoranaPropagationCache}, theta;
    fused::Bool=false, truncate_each_mr=nothing, kwargs...)

    if !fused
        return invoke(PropagationBase.applymergetruncate!,
            Tuple{FermionicRotation,MP.AbstractMajoranaPropagationCache,typeof(theta)},
            gate, prop_cache, theta; truncate_each_mr, kwargs...)
    end

    ms_rotations, coeffs, truncate_after_each_majrot = MP.getmajoranarotations(gate, nsites(prop_cache))
    if !isnothing(truncate_each_mr)
        truncate_after_each_majrot = truncate_each_mr
    end

    for (ii, (gate_ms, coeff)) in enumerate(zip(ms_rotations, coeffs))
        is_last = ii == length(ms_rotations)
        # factor 2 since `MajoranaRotation` implements exp(-i * theta/2 * mstring)
        if truncate_after_each_majrot || is_last
            applymergetruncate!(gate_ms, prop_cache, theta * coeff * 2.0; fused=true, force_truncate=is_last, kwargs...)
        else
            # the rotations in between create strings that only the later rotations pair up again
            applymergetruncate!(gate_ms, prop_cache, theta * coeff * 2.0; fused=true, kwargs..., _NOTRUNC...)
        end
    end

    return prop_cache
end

"""
    applymergetruncate!(gate::ImaginaryFermionicRotation, prop_cache::VectorMajoranaPropagationCache, beta; fused::Bool=false, truncate_each_mr=nothing, normalize_coeffs=true, kwargs...)

Fused overload for `ImaginaryFermionicRotation`, the imaginary-time counterpart of the `FermionicRotation` overload.
With `normalize_coeffs=true` the coefficient truncations are applied after the renormalization, as in the default method, so only the truncations reading the Majorana string are fused into the gate application.
"""
function PropagationBase.applymergetruncate!(gate::ImaginaryFermionicRotation, prop_cache::MP.VectorMajoranaPropagationCache, beta;
    fused::Bool=false, truncate_each_mr=nothing, normalize_coeffs::Bool=true, kwargs...)

    if !fused
        return invoke(PropagationBase.applymergetruncate!,
            Tuple{ImaginaryFermionicRotation,MP.AbstractMajoranaPropagationCache,typeof(beta)},
            gate, prop_cache, beta; truncate_each_mr, normalize_coeffs, kwargs...)
    end

    ms_rotations, coeffs, truncate_after_each_majrot = MP.getmajoranarotations(gate, nsites(prop_cache))
    if !isnothing(truncate_each_mr)
        truncate_after_each_majrot = truncate_each_mr
    end

    for (ii, (gate_ms, coeff)) in enumerate(zip(ms_rotations, coeffs))
        is_last = ii == length(ms_rotations)
        truncate_now = truncate_after_each_majrot || is_last

        if normalize_coeffs
            # the coefficient truncations must see the renormalized coefficients
            applymergetruncate!(gate_ms, prop_cache, beta * coeff; fused=true, truncate_in_merge=false, kwargs..., (truncate_now ? (;) : _NOTRUNC)...)
            mult!(prop_cache, 1 / _identitycoeff(prop_cache))
            truncate_now && truncate!(prop_cache; kwargs...)
        elseif truncate_now
            applymergetruncate!(gate_ms, prop_cache, beta * coeff; fused=true, force_truncate=is_last, kwargs...)
        else
            applymergetruncate!(gate_ms, prop_cache, beta * coeff; fused=true, kwargs..., _NOTRUNC...)
        end
    end

    return prop_cache
end

# Coefficient of the identity string among the active terms. It is the smallest unsigned value, so
# it sits first once the active range is sorted; `getmergedcoeff` would instead scan the inactive
# capacity of the backing arrays when the identity is absent.
function _identitycoeff(prop_cache::MP.VectorMajoranaPropagationCache)
    terms = activeterms(prop_cache)
    coeffs = activecoeffs(prop_cache)
    n = length(terms)
    n == 0 && return zero(eltype(coeffs))
    if sortedprefix(mainsum(prop_cache)) == n
        return iszero(terms[1]) ? coeffs[1] : zero(eltype(coeffs))
    end
    i = findfirst(iszero, terms)
    return isnothing(i) ? zero(eltype(coeffs)) : coeffs[i]
end
