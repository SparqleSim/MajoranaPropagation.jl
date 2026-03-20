function Base.merge!(prop_cache::MajoranaMultiPropagationCache; kwargs...)
    prop_cache = _merge!(prop_cache; kwargs...)
    return prop_cache
end

function _merge!(
    prop_cache::MajoranaMultiPropagationCache{MMS};
    n_levels,
    kwargs...,
) where {TT<:Integer,CT,MMS<:MajoranaSumMulti{TT,CT}}
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    msum_weight_sectors = Set{Int}()
    for key in keys(msum.MultiMajoranas)
        push!(msum_weight_sectors, _get_weight_from_key(key))
    end

    change_weight_sector = [-2, 0, 2]
    for ΔW in change_weight_sector
        for level = 1:n_levels
            for weight_sector in msum_weight_sectors
                dict_key = "$weight_sector-$level"
                ΔW_key = "$(weight_sector + ΔW)-$level"

                if !haskey(aux_msum, dict_key)
                    continue
                end

                aux_dict = aux_msum[dict_key].MultiMajoranas[ΔW_key]
                if isempty(aux_dict)
                    continue
                end
                _merge_single!(msum, aux_dict, ΔW_key)
            end
        end
    end

    for l1 = 1:n_levels
        for l2 = 1:n_levels
            for weight_sector in msum_weight_sectors
                dict_key = "$weight_sector-$l1"
                Δl_key = "$weight_sector-$l2"

                if !haskey(aux_msum, dict_key)
                    continue
                end

                aux_dict = aux_msum[dict_key].MultiMajoranas[Δl_key]
                if isempty(aux_dict)
                    continue
                end

                _merge_single!(msum, aux_dict, Δl_key)
            end
        end
    end

    # remove empty dicts
    for dict_key in collect(keys(msum.MultiMajoranas))
        if isempty(msum.MultiMajoranas[dict_key])
            delete!(msum.MultiMajoranas, dict_key)
        end
    end

    setmainsum!(prop_cache, msum)
    setauxsum!(prop_cache, aux_msum)

    return prop_cache
end

function _merge_single!(msum, aux_dict, dict_key)
    if haskey(msum.MultiMajoranas, dict_key)
        mergewith!(+, msum.MultiMajoranas[dict_key], aux_dict)
        empty!(aux_dict)
    else
        msum.MultiMajoranas[dict_key] = copy(aux_dict)
        empty!(aux_dict)
    end
end

function mergeandempty!(msums::MajoranaSumMulti{TT,CT}, aux_psum; merge_sector=true) where {TT<:Integer,CT}
    prop_cache = MajoranaMultiPropagationCache(msums, aux_psum)
    merge!(prop_cache; merge_sector)
    return mainsum(prop_cache), auxsum(prop_cache)
end

function key_new_weight_sector(key::String, delta_weight::Int)
    weight_sector, level = _split_key(key)
    new_weight_sector = weight_sector + delta_weight
    return "$new_weight_sector-$level"
end