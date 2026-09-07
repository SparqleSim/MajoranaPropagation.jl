#=
Necessary overloads
Is taken 1:1 to the MajoranaPropagation version
=#

abstract type AbstractHybridPropagationCache <: PauliPropagation.AbstractPropagationCache end

mutable struct HybridPropagationCache{HS<:AbstractHybridSum} <: AbstractHybridPropagationCache
    main_hsum::HS
    aux_hsum::HS
end

nfermions(prop_cache::AbstractHybridPropagationCache) = nfermions(mainsum(prop_cache))

HybridPropagationCache(hsum::HS) where {HS<:AbstractHybridSum} = HybridPropagationCache(hsum, similar(hsum))
PauliPropagation.PropagationBase.PropagationCache(hsum::HS) where {HS<:AbstractHybridSum} = HybridPropagationCache(hsum)

terms(prop_cache::HybridPropagationCache) = terms(mainsum(prop_cache))
PauliPropagation.PropagationBase.terms(prop_cache::HybridPropagationCache) = terms(prop_cache)
PauliPropagation.PropagationBase.coefficients(prop_cache::HybridPropagationCache) = coefficients(mainsum(prop_cache))

PauliPropagation.PropagationBase.mainsum(prop_cache::HybridPropagationCache) = prop_cache.main_hsum
PauliPropagation.PropagationBase.auxsum(prop_cache::HybridPropagationCache) = prop_cache.aux_hsum

function PauliPropagation.PropagationBase.setmainsum!(prop_cache::HybridPropagationCache, hsum::HS) where {HS<:AbstractHybridSum}
    prop_cache.main_hsum = hsum
    return prop_cache
end

function PauliPropagation.PropagationBase.setauxsum!(prop_cache::HybridPropagationCache, hsum::HS) where {HS<:AbstractHybridSum}
    prop_cache.aux_hsum = hsum
    return prop_cache
end 