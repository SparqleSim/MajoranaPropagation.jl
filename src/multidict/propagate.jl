function Base.merge!(prop_cache::MajoranaMultiPropagationCache; to, kwargs...)
    prop_cache = _merge_and_empty!(prop_cache; to, kwargs...)
    return prop_cache
end

function _merge_and_empty!(
    prop_cache::MajoranaMultiPropagationCache{MMS};
    level_mapper::Function,
    to,
    kwargs...,
) where {TT<:Integer,CT,MMS<:MajoranaSumMulti{TT,CT}}
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    @timeit to "loop" for iw in eachindex(msum.MultiMajoranas)
        for (ms_int, coeff) in aux_msum.MultiMajoranas[iw]
            add!(msum, level_mapper(ms_int), ms_int, coeff)
        end
        empty!(aux_msum.MultiMajoranas[iw])
    end

    @timeit to "set" setmainsum!(prop_cache, msum)
    @timeit to "set" setauxsum!(prop_cache, aux_msum)

    return prop_cache
end
