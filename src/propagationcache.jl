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

majoranatype(propcache::AbstractMajoranaPropagationCache) = majoranatype(mainsum(propcache))

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
`gate_int`: the tail is `gate_int ⊻ (an ascending subset of the sorted prefix)`, so it is
sorted by `popcount(gate_int)` parallel block-swap passes (see `_xorsorttail!`) instead of
a comparison sort. The two sorted runs are then combined with the same parallel two-pointer
merge kernel as `sortedtailmerge!`. Falls back to the generic `merge!` whenever an XOR
precondition does not hold (no valid sorted prefix, non-CPU storage, or a gate/eltype
mismatch) or when the weight of `gate_int` is above 4.
"""
function xorsortedtailmerge!(prop_cache::VectorMajoranaPropagationCache, gate_int::Union{Integer,Nothing}; thread::Bool=true, truncfunc=nothing, kwargs...)
    n_old = sortedprefix(mainsum(prop_cache))
    n_new = activesize(prop_cache)
    n_tail = n_new - n_old

    if n_tail == 0
        # nothing was appended; the whole active range is still sorted and deduplicated
        return prop_cache
    end

    main_terms, main_coeffs, aux_terms, aux_coeffs = PropagationBase._mainauxarrays(prop_cache)

    if !(main_terms isa Vector{<:Unsigned} && gate_int isa eltype(main_terms)) ||
       n_old <= 0 || n_old > n_new || get_weight(gate_int) > 4 
       @warn "Resorting back to standard `merge!`, n_old=$n_old, n_new=$n_new, w(gate_int)=$get_weight(gate_int)"
        return merge!(prop_cache; thread, truncfunc, kwargs...)
    end

    # ping-pong pair A: the appended tail, in place at the end of the main arrays
    a_terms = view(main_terms, n_old+1:n_new)
    a_coeffs = view(main_coeffs, n_old+1:n_new)

    # ping-pong pair B: as in sortedtailmerge!, the merge output occupies at most aux[1:n_new],
    # so spare aux capacity beyond n_new is free scratch; else allocate fresh. Both branches
    # produce views of Vectors so the A/B swaps in _xorsorttail! stay type-stable.
    if length(aux_terms) - n_new >= n_tail
        b_terms = view(aux_terms, n_new+1:n_new+n_tail)
        b_coeffs = view(aux_coeffs, n_new+1:n_new+n_tail)
    else
        b_terms = view(Base.similar(main_terms, n_tail), 1:n_tail)
        b_coeffs = view(Base.similar(main_coeffs, n_tail), 1:n_tail)
    end

    # popcount(gate_int) parallel block-swap passes; the sorted tail lands back in A when
    # the popcount is even (always, for even-weight Majorana rotation strings), in B when odd
    tail_terms, tail_coeffs = _xorsorttail!(a_terms, a_coeffs, b_terms, b_coeffs, gate_int; thread)

    task_partitioner, n_tasks = PropagationBase._preparetasks(n_old, thread)

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
        write_offsets_per_task = PropagationBase._offsetsfromcounts(merged_counts_per_task)
        merged_count = write_offsets_per_task[end] - 1

        # real pass: each task redoes the same merge, now writing directly into its final position
        AK.itask_partition(n_tasks, n_tasks, 1) do task_id, _
            head_range = task_partitioner[task_id]
            PropagationBase._tailmerge_write!(aux_terms, aux_coeffs, write_offsets_per_task[task_id],
                main_terms, main_coeffs, head_range.start, head_range.stop,
                tail_terms, tail_coeffs, tail_bounds_per_task[task_id], tail_bounds_per_task[task_id+1] - 1, truncfunc, Val(true))
        end
    end

    return PropagationBase._commitwrite!(prop_cache, merged_count, merged_count)
end

"""
    _xorsorttail!(a_terms, a_coeffs, b_terms, b_coeffs, g; thread=true)

Sort the tail held in pair A, given `a_terms[i] == sources[i] ⊻ g` for a strictly
ascending sequence of sources (the anticommuting terms of the sorted prefix, in order),
moving coefficients along with their terms. Pair B is same-length scratch. Returns the
pair holding the sorted tail: A when `popcount(g)` is even (always, for even-weight
Majorana rotation strings), B when odd.

XOR with a fixed `g` preserves the relative order of two values unless the highest bit
in which they differ is a set bit of `g`. So the tail is sorted by one block-swap pass
per set bit `b` of `g`, most significant first: before the pass the buffer is ascending
in `value ⊻ (g masked to bits ≤ b)`, and swapping the two bit-`b` half-blocks within
every run of elements agreeing on the value bits above `b` restores the invariant with
`b` cleared. Each pass is a parallel scatter of contiguous block copies (see
`_xorsortpass!`), so the total cost is `popcount(g)` passes of O(n_tail) instead of a
comparison sort.
"""
function _xorsorttail!(a_terms, a_coeffs, b_terms, b_coeffs, g::TT; thread::Bool=true) where {TT<:Unsigned}
    src_terms, src_coeffs = a_terms, a_coeffs
    dst_terms, dst_coeffs = b_terms, b_coeffs
    remaining_bits = g
    nbits = 8 * sizeof(TT)
    while !iszero(remaining_bits)
        b = (nbits - 1) - leading_zeros(remaining_bits)
        remaining_bits ⊻= one(TT) << b
        _xorsortpass!(dst_terms, dst_coeffs, src_terms, src_coeffs, b; thread)
        src_terms, dst_terms = dst_terms, src_terms
        src_coeffs, dst_coeffs = dst_coeffs, src_coeffs
    end
    return src_terms, src_coeffs
end

# One block-swap pass for set bit b, scattered over contiguous source chunks: each task
# writes the destinations of its own source elements, which are disjoint across tasks
# (sub-blocks map affinely to their destinations and the permutation is a bijection), so
# no synchronization is needed. Separate from the ping-pong loop in _xorsorttail! so the
# closure only captures never-reassigned arguments (no boxing).
function _xorsortpass!(dst_terms, dst_coeffs, src_terms, src_coeffs, b::Int; thread::Bool=true)
    AK.task_partition(length(src_terms), maxtasks(thread), _MIN_ELEMS_PER_TASK) do chunk
        _xorsortpass_chunk!(dst_terms, dst_coeffs, src_terms, src_coeffs, b, first(chunk), last(chunk))
    end
    return
end

function _xorsortpass_chunk!(dst_terms, dst_coeffs, src_terms::AbstractVector{TT}, src_coeffs,
    b::Int, c_lo::Int, c_hi::Int) where {TT<:Unsigned}
    m = length(src_terms)
    hishift = b + 1

    i = c_lo
    @inbounds begin
        # a group straddling the left chunk boundary is recovered in full by binary search
        # (the neighboring task does the same) and only our clipped slice of it is copied
        if c_lo > 1 && (src_terms[c_lo-1] >> hishift) == (src_terms[c_lo] >> hishift)
            i_g = _xorgroupfirst(src_terms, hishift, c_lo)
            j_g = _xorgrouplast(src_terms, hishift, c_lo, m) + 1
            m_split = _xorbitsplit(src_terms, b, i_g, j_g)
            _xormovegroup!(dst_terms, dst_coeffs, src_terms, src_coeffs, i_g, m_split, j_g, c_lo, c_hi)
            i = j_g  # may exceed c_hi when the group swallows the whole chunk
        end

        while i <= c_hi
            # group of elements agreeing on all value bits above b, starting at i:
            # linear scan, counting the leading run with value bit b == 1
            group_hi = src_terms[i] >> hishift
            j = i
            n_first = 0
            while j <= c_hi && (src_terms[j] >> hishift) == group_hi
                n_first += Int(!iszero((src_terms[j] >> b) & one(TT)))
                j += 1
            end
            if j > c_hi && j <= m && (src_terms[j] >> hishift) == group_hi
                # the group continues past the right chunk boundary
                j_g = _xorgrouplast(src_terms, hishift, j, m) + 1
                m_split = _xorbitsplit(src_terms, b, i, j_g)
                _xormovegroup!(dst_terms, dst_coeffs, src_terms, src_coeffs, i, m_split, j_g, c_lo, c_hi)
                break
            end
            # group [i, j) lies fully inside the chunk; its bit-1 elements are the leading n_first
            _xormovegroup!(dst_terms, dst_coeffs, src_terms, src_coeffs, i, i + n_first, j, c_lo, c_hi)
            i = j
        end
    end
    return
end

# block copy with a manual loop below the length where copyto!'s call and bounds-check
# overhead dominates (groups are often singletons, e.g. on low-bit passes of random data)
@inline function _xorcopyrange!(dst_terms, dst_coeffs, src_terms, src_coeffs, d0::Int, s0::Int, n::Int)
    if n < 32
        @inbounds for k in 0:n-1
            dst_terms[d0+k] = src_terms[s0+k]
            dst_coeffs[d0+k] = src_coeffs[s0+k]
        end
    else
        copyto!(dst_terms, d0, src_terms, s0, n)
        copyto!(dst_coeffs, d0, src_coeffs, s0, n)
    end
    return
end

# Move one group's two half-blocks to their swapped destinations -- the bit-1 block
# [i_g, m_split) goes after the bit-0 block [m_split, j_g) -- restricted to the source
# positions in [c_lo, c_hi] owned by the calling task.
@inline function _xormovegroup!(dst_terms, dst_coeffs, src_terms, src_coeffs,
    i_g::Int, m_split::Int, j_g::Int, c_lo::Int, c_hi::Int)
    n_first = m_split - i_g
    n_second = j_g - m_split
    lo = max(i_g, c_lo)
    hi = min(m_split - 1, c_hi)
    if lo <= hi
        _xorcopyrange!(dst_terms, dst_coeffs, src_terms, src_coeffs, lo + n_second, lo, hi - lo + 1)
    end
    lo = max(m_split, c_lo)
    hi = min(j_g - 1, c_hi)
    if lo <= hi
        _xorcopyrange!(dst_terms, dst_coeffs, src_terms, src_coeffs, lo - n_first, lo, hi - lo + 1)
    end
    return
end

# first index in [1, idx] whose value shares the bits above b (i.e. >> hishift) with terms[idx];
# terms >> hishift is non-decreasing, so this is a plain first-true binary search
@inline function _xorgroupfirst(terms::AbstractVector{TT}, hishift::Int, idx::Int) where {TT<:Unsigned}
    @inbounds key = terms[idx] >> hishift
    lo = 1
    hi = idx
    @inbounds while lo < hi
        mid = (lo + hi) >>> 1
        if (terms[mid] >> hishift) == key
            hi = mid
        else
            lo = mid + 1
        end
    end
    return lo
end

# last index in [idx, m] whose value shares the bits above b (i.e. >> hishift) with terms[idx]
@inline function _xorgrouplast(terms::AbstractVector{TT}, hishift::Int, idx::Int, m::Int) where {TT<:Unsigned}
    @inbounds key = terms[idx] >> hishift
    lo = idx
    hi = m
    @inbounds while lo < hi
        mid = (lo + hi + 1) >>> 1
        if (terms[mid] >> hishift) == key
            lo = mid
        else
            hi = mid - 1
        end
    end
    return lo
end

# first index in [i_g, j_g) whose value has bit b == 0 (within a group the bit-1 run comes
# first), or j_g when every element has bit b == 1
@inline function _xorbitsplit(terms::AbstractVector{TT}, b::Int, i_g::Int, j_g::Int) where {TT<:Unsigned}
    lo = i_g
    hi = j_g
    @inbounds while lo < hi
        mid = (lo + hi) >>> 1
        if iszero((terms[mid] >> b) & one(TT))
            hi = mid
        else
            lo = mid + 1
        end
    end
    return lo
end
