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