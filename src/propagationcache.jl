abstract type AbstractMajoranaPropagationCache <: AbstractPropagationCache end

mutable struct MajoranaPropagationCache{MS<:AbstractMajoranaSum} <: AbstractMajoranaPropagationCache
    mainsum::MS
    auxsum::MS
end

MajoranaPropagationCache(msum::MS) where {MS<:AbstractMajoranaSum} = MajoranaPropagationCache(msum, similar(msum))
PropagationBase.PropagationCache(msum::MS) where {MS<:AbstractMajoranaSum} = MajoranaPropagationCache(msum)

PropagationBase.nsites(prop_cache::MajoranaPropagationCache) = nsites(mainsum(prop_cache))
majoranas(prop_cache::MajoranaPropagationCache) = majoranas(mainsum(prop_cache))
PropagationBase.terms(prop_cache::MajoranaPropagationCache) = majoranas(prop_cache)
PropagationBase.coefficients(prop_cache::MajoranaPropagationCache) = coefficients(mainsum(prop_cache))

PropagationBase.mainsum(prop_cache::MajoranaPropagationCache) = prop_cache.mainsum
PropagationBase.auxsum(prop_cache::MajoranaPropagationCache) = prop_cache.auxsum

function PropagationBase.setmainsum!(prop_cache::MajoranaPropagationCache, msum::MS) where {MS<:AbstractMajoranaSum}
    prop_cache.mainsum = msum
    return prop_cache
end

function PropagationBase.setauxsum!(prop_cache::MajoranaPropagationCache, msum::MS) where {MS<:AbstractMajoranaSum}
    prop_cache.auxsum = msum
    return prop_cache
end