abstract type AbstractMajoranaPropagationCache <: AbstractPropagationCache end

nfermions(prop_cache::AbstractMajoranaPropagationCache) = nfermions(mainsum(prop_cache))


mutable struct MajoranaPropagationCache{MS<:AbstractMajoranaSum} <: AbstractMajoranaPropagationCache
    main_msum::MS
    aux_msum::MS
end

MajoranaPropagationCache(msum::MS) where {MS<:AbstractMajoranaSum} = MajoranaPropagationCache(msum, similar(msum))
PropagationBase.PropagationCache(msum::MS) where {MS<:AbstractMajoranaSum} = MajoranaPropagationCache(msum)

majoranas(prop_cache::MajoranaPropagationCache) = majoranas(mainsum(prop_cache))
PropagationBase.terms(prop_cache::MajoranaPropagationCache) = majoranas(prop_cache)
PropagationBase.coefficients(prop_cache::MajoranaPropagationCache) = coefficients(mainsum(prop_cache))

PropagationBase.mainsum(prop_cache::MajoranaPropagationCache) = prop_cache.main_msum
PropagationBase.auxsum(prop_cache::MajoranaPropagationCache) = prop_cache.aux_msum

function PropagationBase.setmainsum!(prop_cache::AbstractMajoranaPropagationCache, msum::MS) where {MS<:AbstractMajoranaSum}
    prop_cache.main_msum = msum
    return prop_cache
end

function PropagationBase.setauxsum!(prop_cache::AbstractMajoranaPropagationCache, msum::MS) where {MS<:AbstractMajoranaSum}
    prop_cache.aux_msum = msum
    return prop_cache
end

# VectorMajoranaPropagationCache
mutable struct VectorMajoranaPropagationCache{VMS<:VectorMajoranaSum,VB,VI} <: AbstractMajoranaPropagationCache
    main_msum::VMS
    aux_msum::VMS
    flags::VB
    indices::VI
    active_size::Int
end

# Overload for generality
function PropagationBase.PropagationCache(vecmsum::VectorMajoranaSum)
    return VectorMajoranaPropagationCache(vecmsum)
end

function VectorMajoranaPropagationCache(vecmsum::VectorMajoranaSum{VT,VC}) where {VT,VC}
    aux_vecmsum = Base.similar(vecmsum)
    flags = Base.similar(majoranas(vecmsum), Bool)
    indices = Base.similar(majoranas(vecmsum), Int)
    return VectorMajoranaPropagationCache(vecmsum, aux_vecmsum, flags, indices, length(vecmsum))
end

PropagationBase.mainsum(vprop_cache::VectorMajoranaPropagationCache) = vprop_cache.main_msum
PropagationBase.auxsum(vprop_cache::VectorMajoranaPropagationCache) = vprop_cache.aux_msum

function VectorMajoranaPropagationCache(msum::MajoranaSum)
    return VectorMajoranaPropagationCache(VectorMajoranaSum(msum))
end

# Convert back to vector and dense sums
function VectorMajoranaSum(prop_cache::VectorMajoranaPropagationCache)
    vecmsum = deepcopy(mainsum(prop_cache))
    resize!(vecmsum, activesize(prop_cache))
    return vecmsum
end

function MajoranaSum(prop_cache::VectorMajoranaPropagationCache)
    merge!(prop_cache)
    return MajoranaSum(nqubits(prop_cache), Dict(zip(activeterms(prop_cache), activecoeffs(prop_cache))))
end

PropagationBase.activesize(prop_cache::VectorMajoranaPropagationCache) = prop_cache.active_size
PropagationBase.setactivesize!(prop_cache::VectorMajoranaPropagationCache, new_size::Int) = (prop_cache.active_size = new_size; prop_cache)

PropagationBase.indices(prop_cache::VectorMajoranaPropagationCache) = prop_cache.indices
PropagationBase.flags(prop_cache::VectorMajoranaPropagationCache) = prop_cache.flags

# Term and coefficient accessors for caches
majoranas(prop_cache::VectorMajoranaPropagationCache) = activeterms(prop_cache)
PropagationBase.coefficients(prop_cache::VectorMajoranaPropagationCache) = activecoeffs(prop_cache)

function Base.resize!(prop_cache::VectorMajoranaPropagationCache, n_new::Int)
    resize!(prop_cache.main_msum, n_new)
    resize!(prop_cache.aux_msum, n_new)
    resize!(prop_cache.flags, n_new)
    resize!(prop_cache.indices, n_new)
    return prop_cache
end

# ========== sorted-runs merge ========== #

# Called once at the start of propagate! so that the active terms are sorted and
# deduplicated. From then on, every merge only has to deal with the unsorted tail
# appended by applytoall! (see mergesortedruns!).
_presortcache!(prop_cache::AbstractMajoranaPropagationCache) = prop_cache
_presortcache!(prop_cache::VectorMajoranaPropagationCache) = merge!(prop_cache)

# The size of the sorted, deduplicated prefix before applytoall! appends new terms.
# For dict caches there is no such notion; merging is cheap there anyway.
_sortedprefixsize(prop_cache::AbstractMajoranaPropagationCache) = 0
_sortedprefixsize(prop_cache::VectorMajoranaPropagationCache) = activesize(prop_cache)

_mergeafterapply!(prop_cache::AbstractMajoranaPropagationCache, n_sorted::Int, gate_int; kwargs...) = merge!(prop_cache; kwargs...)
_mergeafterapply!(prop_cache::VectorMajoranaPropagationCache, n_sorted::Int, gate_int; kwargs...) = mergesortedruns!(prop_cache, n_sorted, gate_int; kwargs...)

"""
    mergesortedruns!(prop_cache::VectorMajoranaPropagationCache, n_sorted::Int, [gate_int]; kwargs...)

Merge and deduplicate the active terms of `prop_cache`, exploiting that the first `n_sorted`
active terms are already sorted and deduplicated (the state a previous `merge!` or
`mergesortedruns!` leaves behind, which `applytoall!` and `truncate!` preserve), while the
remaining active terms are the unsorted tail appended by `applytoall!`.

Instead of the full `sortperm` over all active terms that `merge!` performs, this only sorts
the tail and then combines the two sorted runs in a single linear pass, merging coefficients
of equal terms on the fly. Falls back to a full `merge!` if the sorted-prefix assumption
cannot be verified.

If the Majorana string `gate_int` of the rotation that produced the tail is passed, the tail
permutation is computed with `popcount(gate_int)` linear block-swap passes (see
`_xorsortperm!`) instead of a comparison sort.
"""
mergesortedruns!(prop_cache::VectorMajoranaPropagationCache, n_sorted::Int; kwargs...) =
    mergesortedruns!(prop_cache, n_sorted, nothing; kwargs...)

function mergesortedruns!(prop_cache::VectorMajoranaPropagationCache, n_sorted::Int, gate_int::Union{Integer,Nothing}; kwargs...)
    n_total = activesize(prop_cache)

    if n_total == 0
        return prop_cache
    end

    main_terms = majoranas(mainsum(prop_cache))
    main_coeffs = coefficients(mainsum(prop_cache))

    # the fast path requires CPU arrays (scalar indexing) and a valid sorted prefix
    if !(main_terms isa Vector) ||
       n_sorted <= 0 ||
       n_sorted > n_total ||
       !issorted(view(main_terms, 1:n_sorted))
        return merge!(prop_cache; kwargs...)
    end

    m = n_total - n_sorted
    if m == 0
        # fully sorted and already deduplicated by the previous merge
        return prop_cache
    end

    aux_terms = majoranas(auxsum(prop_cache))
    aux_coeffs = coefficients(auxsum(prop_cache))
    @assert length(aux_terms) >= n_total "VectorMajoranaPropagationCache aux terms array is not large enough to hold the merged terms."

    tail_terms = view(main_terms, n_sorted+1:n_total)
    tail_coeffs = view(main_coeffs, n_sorted+1:n_total)

    # sort a permutation of the tail, reusing the cache's indices array as scratch
    tail_perm = view(prop_cache.indices, 1:m)
    if gate_int isa eltype(main_terms) &&
       eltype(main_terms) <: Unsigned &&
       m <= _XORSORT_MAX_TAIL &&
       length(prop_cache.indices) >= 2 * m
        # the tail is gate_int ⊻ (terms in ascending order), so its sorted order follows
        # from popcount(gate_int) linear block-swap passes instead of a comparison sort
        scratch = view(prop_cache.indices, m+1:2*m)
        _xorsortperm!(tail_perm, scratch, tail_terms, gate_int)
    else
        AK.sortperm!(tail_perm, tail_terms)
    end

    # single linear pass over both sorted runs, writing into the aux arrays;
    # equal terms end up adjacent in the output, so coefficients are combined on the fly
    k = 0
    i = 1
    j = 1
    @inbounds while i <= n_sorted || j <= m
        take_prefix = j > m || (i <= n_sorted && !isless(tail_terms[tail_perm[j]], main_terms[i]))

        if take_prefix
            term = main_terms[i]
            coeff = main_coeffs[i]
            i += 1
        else
            p = tail_perm[j]
            term = tail_terms[p]
            coeff = tail_coeffs[p]
            j += 1
        end

        if k > 0 && aux_terms[k] == term
            aux_coeffs[k] = mergefunc(aux_coeffs[k], coeff)
        else
            k += 1
            aux_terms[k] = term
            aux_coeffs[k] = coeff
        end
    end

    swapsums!(prop_cache)
    setactivesize!(prop_cache, k)

    return prop_cache
end

# Above this tail size, the multithreaded comparison sort outpaces the serial
# XOR block-swap shuffle (both are memory-bandwidth-bound at large sizes, but the
# sample sort uses all threads).
const _XORSORT_MAX_TAIL = 100_000

"""
    _xorsortperm!(perm, scratch, tail_terms, g)

Fill `perm` with the permutation that sorts `tail_terms`, given that
`tail_terms[i] == sources[i] ⊻ g` for a strictly ascending sequence of sources
(the anticommuting terms of the sorted prefix, in order).

XOR with a fixed `g` preserves the relative order of two values unless the highest bit
in which they differ is a set bit of `g`. Consequently the sorted order of the tail is
obtained by, for each set bit `b` of `g` from most significant to least significant,
swapping the two half-blocks (bit `b` = 0 and bit `b` = 1 of the source) within every
group of elements that agree on all source bits above `b`. Each pass is a single linear
gather, so the total cost is `popcount(g)` passes of O(length(perm)) instead of a
comparison sort. `scratch` must have the same length as `perm`.
"""
function _xorsortperm!(perm::AbstractVector{Int}, scratch::AbstractVector{Int}, tail_terms, g::TT) where {TT<:Integer}
    m = length(perm)

    @inbounds for i in 1:m
        perm[i] = i
    end

    src = perm
    dst = scratch
    remaining_bits = g
    nbits = 8 * sizeof(TT)

    @inbounds while !iszero(remaining_bits)
        b = (nbits - 1) - leading_zeros(remaining_bits)
        remaining_bits ⊻= one(TT) << b
        hishift = b + 1

        i = 1
        while i <= m
            # scan the group of elements agreeing on all source bits above b,
            # counting how many have source bit b == 0 (they come first, sources ascend)
            group_hi = (tail_terms[src[i]] ⊻ g) >> hishift
            j = i
            n0 = 0
            while j <= m
                source = tail_terms[src[j]] ⊻ g
                (source >> hishift) == group_hi || break
                n0 += Int(iszero((source >> b) & one(TT)))
                j += 1
            end

            # bit b of the image is flipped, so the bit-1 sources come first in the output
            n1 = (j - i) - n0
            copyto!(dst, i, src, i + n0, n1)
            copyto!(dst, i + n1, src, i, n0)
            i = j
        end

        src, dst = dst, src
    end

    if src !== perm
        copyto!(perm, src)
    end
    return perm
end