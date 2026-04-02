# multi dict propagation cache

mutable struct MajoranaMultiPropagationCache{MMS<:MajoranaSumMulti} <: AbstractMajoranaPropagationCache
    main_msum::MMS
    aux_msum::MMS
end

# Overload for generality
function PropagationBase.PropagationCache(multimsum::MajoranaSumMulti)
    return MajoranaMultiPropagationCache(multimsum)
end

function MajoranaMultiPropagationCache(multimsum::MajoranaSumMulti{TT,VC}) where {TT<:Integer,VC}
    aux_msum = similar(multimsum)
    return MajoranaMultiPropagationCache(multimsum, aux_msum)
end

PropagationBase.mainsum(multimsum::MajoranaMultiPropagationCache) = multimsum.main_msum
PropagationBase.auxsum(multimsum::MajoranaMultiPropagationCache) = multimsum.aux_msum
nfermions(prop_cache::MajoranaMultiPropagationCache) = nfermions(mainsum(prop_cache))

function PropagationBase.setmainsum!(prop_cache::MajoranaMultiPropagationCache, msum::MajoranaSumMulti)
    prop_cache.main_msum = msum
    return prop_cache
end

function PropagationBase.setauxsum!(
    prop_cache::MajoranaMultiPropagationCache,
    aux_msum::MajoranaSumMulti,
)
    prop_cache.aux_msum = aux_msum
    return prop_cache
end

function MajoranaMultiPropagationCache(msum::MajoranaSum)
    return MajoranaMultiPropagationCache(MajoranaSumMulti(msum))
end

function show_stats(prop_cache::MajoranaMultiPropagationCache{MMS}) where {MMS<:MajoranaSumMulti}
    #println("Main sum stats:")
    show_stats(mainsum(prop_cache))
    #println("Auxiliary sum stats:")
    #show_stats(auxsum(prop_cache))
end

Base.length(prop_cache::MajoranaMultiPropagationCache) = length(mainsum(prop_cache))