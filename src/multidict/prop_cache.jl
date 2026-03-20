# multi dict propagation cache

mutable struct MajoranaMultiPropagationCache{MMS<:MajoranaSumMulti} <: AbstractMajoranaPropagationCache
    main_msum::MMS
    aux_msum::Dict{Int64,MMS}
end

# Overload for generality
function PropagationBase.PropagationCache(multimsum::MajoranaSumMulti)
    return MajoranaMultiPropagationCache(multimsum)
end

function MajoranaMultiPropagationCache(multimsum::MajoranaSumMulti{TT,VC}, nlevels::Int) where {TT<:Integer,VC}
    all_aux_msums::Dict{Int64,MajoranaSumMulti{TT,VC}} = Dict()
    #@show multimsum
    for k in keys(multimsum)
        weight_key = _get_weight_from_key(k, nlevels)
        all_aux_msums[k] = similar(multimsum, weight_key, nlevels)
    end
    return MajoranaMultiPropagationCache(multimsum, all_aux_msums)
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
    aux_msum::Dict{Int64,<:MajoranaSumMulti},
)
    prop_cache.aux_msum = aux_msum
    return prop_cache
end

function MajoranaMultiPropagationCache(msum::MajoranaSum)
    return MajoranaMultiPropagationCache(MajoranaSumMulti(msum))
end
