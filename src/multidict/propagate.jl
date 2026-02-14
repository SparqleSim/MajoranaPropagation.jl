function Base.merge!(prop_cache::MajoranaMultiPropagationCache; merge_sector=true, kwargs...)
    prop_cache = _merge!(prop_cache; merge_sector, kwargs...)
    return prop_cache
end

function _merge!(
    prop_cache::MajoranaMultiPropagationCache{MMS};
    merge_sector=true,
    kwargs...,
) where {TT<:Integer,CT,MMS<:MajoranaSumMulti{TT,CT}}
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    sorted_keys = sort(collect(keys(msum)))

    # do the Delta W = 0 merges first
    if merge_sector
        @threads for weight_key in sorted_keys
            if haskey(aux_msum, weight_key)
                mergewith!(+, msum.MultiMajoranas[weight_key], aux_msum[weight_key].MultiMajoranas[weight_key])
                empty!(aux_msum[weight_key].MultiMajoranas[weight_key])
            end
        end
    end

    # do the Delta W = +2 merges
    @threads for weight_key in sorted_keys
        if haskey(aux_msum, weight_key)
            aux_dict = aux_msum[weight_key].MultiMajoranas
            if !isempty(aux_dict[weight_key + 2])
                if !haskey(msum.MultiMajoranas, weight_key + 2)
                    msum.MultiMajoranas[weight_key + 2] = Dict{TT,CT}()
                end
                mergewith!(+, msum.MultiMajoranas[weight_key + 2], aux_dict[weight_key + 2])
                empty!(aux_dict[weight_key + 2])
            end
        end
    end

    # do the Delta W = -2 merges
    @threads for weight_key in sorted_keys
        if haskey(aux_msum, weight_key)
            aux_dict = aux_msum[weight_key].MultiMajoranas
            if !isempty(aux_dict[weight_key - 2])
                if !haskey(msum.MultiMajoranas, weight_key - 2)
                    msum.MultiMajoranas[weight_key - 2] = Dict{TT,CT}()
                end
                mergewith!(+, msum.MultiMajoranas[weight_key - 2], aux_dict[weight_key - 2])
                empty!(aux_dict[weight_key - 2])
            end
        end
    end

    # remove empty dicts
    for weight_key in sort(collect(keys(msum.MultiMajoranas)))
        if isempty(msum.MultiMajoranas[weight_key])
            delete!(msum.MultiMajoranas, weight_key)
        end
    end

    setmainsum!(prop_cache, msum)
    setauxsum!(prop_cache, aux_msum)

    return prop_cache
end

function mergeandempty!(msums::MajoranaSumMulti{TT,CT}, aux_psum; merge_sector=true) where {TT<:Integer,CT}
    prop_cache = MajoranaMultiPropagationCache(msums, aux_psum)
    merge!(prop_cache; merge_sector)
    return mainsum(prop_cache), auxsum(prop_cache)
end