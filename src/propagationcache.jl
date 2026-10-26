abstract type AbstractMajoranaPropagationCache <: AbstractPropagationCache end

"""
    nfermions(prop_cache::AbstractMajoranaPropagationCache)

Get the number of fermions that the main Majorana sum of the propagation cache is defined on.
"""
nfermions(prop_cache::AbstractMajoranaPropagationCache) = nfermions(mainsum(prop_cache))


"""
    MajoranaPropagationCache(msum::MS) where {MS<:AbstractMajoranaSum}

A cache structure used during propagation, holding the main Majorana sum `main_msum` and an auxiliary sum `aux_msum` into which newly created Majorana strings are written before being merged.
"""
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

    # the array kernels index all four arrays up to the active size without bounds checks, and
    # write at positions they accumulate in `indices`, which must not wrap
    function VectorMajoranaPropagationCache(main_msum::VMS, aux_msum::VMS, flags::VB, indices::VI, active_size::Int) where {VMS<:VectorMajoranaSum,VB,VI}
        n = length(main_msum)
        if !(length(aux_msum) == length(flags) == length(indices) == n)
            throw(ArgumentError("main_msum, aux_msum, flags and indices must have the same length."))
        end
        if eltype(indices) !== Int
            throw(ArgumentError("indices must hold Int, got $(eltype(indices))."))
        end
        if !(0 <= active_size <= n)
            throw(ArgumentError("active_size must be between 0 and the length of main_msum, got $active_size for $n."))
        end
        return new{VMS,VB,VI}(main_msum, aux_msum, flags, indices, active_size)
    end
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
"""
    VectorMajoranaSum(prop_cache::VectorMajoranaPropagationCache)

Extract a copy of the active part of the main `VectorMajoranaSum` from the propagation cache.
"""
function VectorMajoranaSum(prop_cache::VectorMajoranaPropagationCache)
    vecmsum = deepcopy(mainsum(prop_cache))
    resize!(vecmsum, activesize(prop_cache))
    return vecmsum
end

"""
    MajoranaSum(prop_cache::VectorMajoranaPropagationCache; thread=true)

Merge the propagation cache and convert its main sum into a dictionary-based `MajoranaSum`.
"""
function MajoranaSum(prop_cache::VectorMajoranaPropagationCache; thread::Bool=true)
    merge!(prop_cache; thread)
    vmsum = mainsum(prop_cache)
    return MajoranaSum(vmsum.nsites, vmsum.is_spinful, Dict(zip(activeterms(prop_cache), activecoeffs(prop_cache))))
end

PropagationBase.activesize(prop_cache::VectorMajoranaPropagationCache) = prop_cache.active_size
function PropagationBase.setactivesize!(prop_cache::VectorMajoranaPropagationCache, new_size::Int)
    if !(0 <= new_size <= capacity(prop_cache))
        throw(ArgumentError("active size must be between 0 and the capacity, got $new_size for $(capacity(prop_cache))."))
    end
    prop_cache.active_size = new_size
    return prop_cache
end

PropagationBase.indices(prop_cache::VectorMajoranaPropagationCache) = prop_cache.indices
PropagationBase.flags(prop_cache::VectorMajoranaPropagationCache) = prop_cache.flags

# the active terms as a sum of views, without the room past them
function PropagationBase.activesum(prop_cache::VectorMajoranaPropagationCache)
    vmsum = mainsum(prop_cache)
    n_sorted = sortedprefix(vmsum)
    # a sorted prefix that reaches past the active terms cannot be trusted
    n_sorted = n_sorted > activesize(prop_cache) ? 0 : n_sorted
    return VectorMajoranaSum(nsites(vmsum), is_spinful(vmsum), activeterms(prop_cache), activecoeffs(prop_cache), n_sorted)
end

# Term and coefficient accessors for caches
majoranas(prop_cache::VectorMajoranaPropagationCache) = activeterms(prop_cache)
PropagationBase.coefficients(prop_cache::VectorMajoranaPropagationCache) = activecoeffs(prop_cache)

"""
    resize!(prop_cache::VectorMajoranaPropagationCache, n_new::Int)

Resize all internal arrays of the propagation cache to length `n_new`, keeping at most `n_new` terms active.
"""
function Base.resize!(prop_cache::VectorMajoranaPropagationCache, n_new::Int)
    # the sums grow their own arrays, and the cache its flags and indices
    resize!(prop_cache.main_msum, n_new)
    resize!(prop_cache.aux_msum, n_new)
    prop_cache.flags = PropagationBase._resizearray(prop_cache.flags, n_new)
    prop_cache.indices = PropagationBase._resizearray(prop_cache.indices, n_new)
    setactivesize!(prop_cache, min(activesize(prop_cache), n_new))
    return prop_cache
end
