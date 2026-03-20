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
        push!(msum_weight_sectors, _get_weight_from_key(key, n_levels))
    end

    change_weight_sector = [-2, 0, 2]
    for ΔW in change_weight_sector
        weight_sectors_vec = collect(msum_weight_sectors)
        # Thread-local storage for accumulated results
        thread_results = [Dict{Int64,Dict{TT,CT}}() for _ in 1:Threads.maxthreadid()]
        
        Threads.@threads for level = 1:n_levels
            tid = Threads.threadid()
            for weight_sector in weight_sectors_vec
                dict_key = Int64((weight_sector - 1) * n_levels + level)
                ΔW_key = Int64((weight_sector + ΔW - 1) * n_levels + level)

                if !haskey(aux_msum, dict_key)
                    continue
                end

                aux_dict = aux_msum[dict_key].MultiMajoranas[ΔW_key]
                if isempty(aux_dict)
                    continue
                end
                # Stage by reference; actual merge + clear happens serially below.
                _merge_single_to_dict!(thread_results[tid], aux_dict, ΔW_key)
            end
        end
        
        # Merge thread results back into msum (single-threaded)
        for tid in eachindex(thread_results)
            for (key, val) in thread_results[tid]
                dest = get!(msum.MultiMajoranas, key, Dict{TT,CT}())
                mergewith!(+, dest, val)
                empty!(val)
            end
        end
    end

    for l1 = 1:n_levels
        weight_sectors_vec = collect(msum_weight_sectors)
        thread_results = [Dict{Int64,Dict{TT,CT}}() for _ in 1:Threads.maxthreadid()]
        
        Threads.@threads for l2 = 1:n_levels
            tid = Threads.threadid()
            for weight_sector in weight_sectors_vec
                dict_key = Int64((weight_sector - 1) * n_levels + l1)
                Δl_key = Int64((weight_sector - 1) * n_levels + l2)

                if !haskey(aux_msum, dict_key)
                    continue
                end

                aux_dict = aux_msum[dict_key].MultiMajoranas[Δl_key]
                if isempty(aux_dict)
                    continue
                end

                _merge_single_to_dict!(thread_results[tid], aux_dict, Δl_key)
            end
        end
        
        for tid in eachindex(thread_results)
            for (key, val) in thread_results[tid]
                dest = get!(msum.MultiMajoranas, key, Dict{TT,CT}())
                mergewith!(+, dest, val)
                empty!(val)
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

function _merge_single_to_dict!(target_dict::Dict{Int64,Dict{TT,CT}}, aux_dict::Dict{TT,CT}, dict_key::Int64) where {TT<:Integer,CT}
    if haskey(target_dict, dict_key)
        mergewith!(+, target_dict[dict_key], aux_dict)
    else
        target_dict[dict_key] = aux_dict
    end
end

function mergeandempty!(msums::MajoranaSumMulti{TT,CT}, aux_psum; merge_sector=true) where {TT<:Integer,CT}
    prop_cache = MajoranaMultiPropagationCache(msums, aux_psum)
    merge!(prop_cache; merge_sector)
    return mainsum(prop_cache), auxsum(prop_cache)
end

function key_new_weight_sector(key::Int64, delta_weight::Int, n_levels::Int)
    weight_sector, level = _split_key(key, n_levels)
    new_weight_sector = weight_sector + delta_weight
    return ((new_weight_sector - 1) * n_levels + level)
end