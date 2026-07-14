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

# ========== gate-aware sorted-tail merge ========== #

# Merge the tail appended by a Majorana rotation into the sorted prefix tracked on the
# vector sum (see `sortedprefix`). Dict caches just merge; vector caches use the gate's
# Majorana string to sort the tail without a comparison sort where possible, and defer to
# the generic `merge!` (which picks between `sortedtailmerge!` and a full sort) otherwise.
_mergeafterapply!(prop_cache::AbstractMajoranaPropagationCache, gate_int; kwargs...) = merge!(prop_cache; kwargs...)
_mergeafterapply!(prop_cache::VectorMajoranaPropagationCache, gate_int; kwargs...) = xorsortedtailmerge!(prop_cache, gate_int; kwargs...)

"""
    xorsortedtailmerge!(prop_cache::VectorMajoranaPropagationCache, gate_int; thread=true, kwargs...)

`PropagationBase.sortedtailmerge!` specialized to a tail appended by the Majorana rotation
`gate_int`: the tail is `gate_int ⊻ (an ascending subset of the sorted prefix)`, so its
sorted permutation follows from `popcount(gate_int)` linear block-swap passes (see
`_xorsortperm!`) instead of a comparison sort. The two sorted runs are then combined with
the same parallel two-pointer merge kernel as `sortedtailmerge!`. Falls back to the generic
`merge!` whenever an XOR precondition does not hold (no valid sorted prefix, non-CPU
storage, tail larger than `_XORSORT_MAX_TAIL`, or insufficient index scratch).
"""
function xorsortedtailmerge!(prop_cache::VectorMajoranaPropagationCache, gate_int::Union{Integer,Nothing}; thread::Bool=true, truncfunc=nothing, kwargs...)
    main_terms = majoranas(mainsum(prop_cache))
    main_coeffs = coefficients(mainsum(prop_cache))
    n_old = sortedprefix(mainsum(prop_cache))
    n_new = activesize(prop_cache)
    n_tail = n_new - n_old

    if n_new == 0
        return prop_cache
    end

    if !(main_terms isa Vector) ||
       !(gate_int isa eltype(main_terms)) ||
       !(eltype(main_terms) <: Unsigned) ||
       n_old <= 0 ||
       n_old > n_new ||
       n_tail > _XORSORT_MAX_TAIL ||
       length(indices(prop_cache)) < 2 * n_tail
        return merge!(prop_cache; thread, truncfunc, kwargs...)
    end

    if n_tail == 0
        # nothing was appended; the whole active range is still sorted and deduplicated
        return prop_cache
    end

    aux_terms = majoranas(auxsum(prop_cache))
    aux_coeffs = coefficients(auxsum(prop_cache))

    unsorted_tail_terms = view(main_terms, n_old+1:n_new)
    unsorted_tail_coeffs = view(main_coeffs, n_old+1:n_new)

    # XOR block-swap permutation of the tail, reusing the cache's indices array as perm + scratch
    tail_perm = view(indices(prop_cache), 1:n_tail)
    perm_scratch = view(indices(prop_cache), n_tail+1:2*n_tail)
    _xorsortperm!(tail_perm, perm_scratch, unsorted_tail_terms, gate_int)

    # as in sortedtailmerge!: the merge output occupies at most aux[1:n_new], so spare aux
    # capacity beyond n_new is free scratch for the sorted tail; else allocate fresh
    if length(aux_terms) - n_new >= n_tail
        tail_terms = view(aux_terms, n_new+1:n_new+n_tail)
        tail_coeffs = view(aux_coeffs, n_new+1:n_new+n_tail)
    else
        tail_terms = Base.similar(unsorted_tail_terms)
        tail_coeffs = Base.similar(unsorted_tail_coeffs)
    end
    permuteviaindices!(tail_terms, tail_coeffs, unsorted_tail_terms, unsorted_tail_coeffs, tail_perm; thread)

    task_partitioner = AK.TaskPartitioner(n_old, maxtasks(thread), _MIN_ELEMS_PER_TASK)
    n_tasks = task_partitioner.num_tasks

    if n_tasks == 1
        merged_count = PropagationBase._tailmerge_write!(aux_terms, aux_coeffs, 1,
            main_terms, main_coeffs, 1, n_old, tail_terms, tail_coeffs, 1, n_tail, truncfunc, Val(true))
    else
        # slice and partition the two-pointer merge across threads (same scheme as sortedtailmerge!)
        tail_bounds_per_task = Vector{Int}(undef, n_tasks + 1)
        tail_bounds_per_task[1] = 1
        tail_bounds_per_task[n_tasks+1] = n_tail + 1
        @inbounds for task_id in 1:(n_tasks-1)
            head_chunk_boundary_term = main_terms[task_partitioner[task_id].stop]
            tail_bounds_per_task[task_id+1] = searchsortedlast(tail_terms, head_chunk_boundary_term) + 1
        end

        # dry run: each task counts its own merged output size (unknown ahead of time due to collisions)
        merged_counts_per_task = Vector{Int}(undef, n_tasks)
        AK.itask_partition(n_tasks, n_tasks, 1) do task_id, _
            head_range = task_partitioner[task_id]
            merged_counts_per_task[task_id] = PropagationBase._tailmerge_write!(aux_terms, aux_coeffs, 1,
                main_terms, main_coeffs, head_range.start, head_range.stop,
                tail_terms, tail_coeffs, tail_bounds_per_task[task_id], tail_bounds_per_task[task_id+1] - 1, truncfunc, Val(false))
        end

        # prefix sum over the per-task counts gives each task its exact final write offset
        write_offsets_per_task = Vector{Int}(undef, n_tasks + 1)
        write_offsets_per_task[1] = 1
        @inbounds for task_id in 1:n_tasks
            write_offsets_per_task[task_id+1] = write_offsets_per_task[task_id] + merged_counts_per_task[task_id]
        end
        merged_count = write_offsets_per_task[n_tasks+1] - 1

        # real pass: each task redoes the same merge, now writing directly into its final position
        AK.itask_partition(n_tasks, n_tasks, 1) do task_id, _
            head_range = task_partitioner[task_id]
            PropagationBase._tailmerge_write!(aux_terms, aux_coeffs, write_offsets_per_task[task_id],
                main_terms, main_coeffs, head_range.start, head_range.stop,
                tail_terms, tail_coeffs, tail_bounds_per_task[task_id], tail_bounds_per_task[task_id+1] - 1, truncfunc, Val(true))
        end
    end

    swapsums!(prop_cache)
    setactivesize!(prop_cache, merged_count)
    setsortedprefix!(mainsum(prop_cache), merged_count)

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