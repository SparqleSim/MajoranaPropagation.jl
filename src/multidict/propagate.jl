function Base.merge!(prop_cache::MajoranaMultiPropagationCache; to=TimerOutput(), kwargs...)
    prop_cache = _merge_and_empty!(prop_cache; to, kwargs...)
    return prop_cache
end

function _merge_and_empty!(
    prop_cache::MajoranaMultiPropagationCache{MMS};
    level_mapper=nothing,
    kwargs...,
) where {TT<:Integer,CT,MMS<:MajoranaSumMulti{TT,CT}}
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)  # Vector{MMS}
    main_levels = msum.MultiMajoranas
    zero_coeff = zero(CT)
    n_levels = length(main_levels)
    _ = level_mapper

    @threads for level in 1:n_levels
        dest_dict = main_levels[level]
        for iw in eachindex(aux_msum)
            src_dict = aux_msum[iw].MultiMajoranas[level]
            isempty(src_dict) && continue

            for (ms_int, coeff) in src_dict
                dest_dict[ms_int] = get(dest_dict, ms_int, zero_coeff) + coeff
            end
            empty!(src_dict)
        end
    end

    setmainsum!(prop_cache, msum)
    setauxsum!(prop_cache, aux_msum)
    return prop_cache
end
