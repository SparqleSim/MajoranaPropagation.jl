function PropagationBase.truncate!(
    prop_cache::MajoranaMultiPropagationCache;
    max_weight::Real=Inf, min_abs_coeff=1e-10, max_unpaired::Real=Inf,
    max_freq::Real=Inf, max_sins::Real=Inf,
    unpaired_mask=nothing,
    customtruncfunc=nothing,
    kwargs...
)
    if isnothing(unpaired_mask)
        unpaired_mask = create_unpaired_mask(nfermions(mainsum(prop_cache)))
    end
    function truncfunc(mstr, coeff)
        is_truncated = false
        if PauliPropagation.truncatemincoeff(coeff, min_abs_coeff)
            is_truncated = true
        elseif truncateunpaired(mstr, max_unpaired, unpaired_mask)
            is_truncated = true
        elseif truncatemajoranaweight(mstr, max_weight)
            is_truncated = true
        elseif PauliPropagation.truncatefrequency(coeff, max_freq)
            is_truncated = true
        elseif PauliPropagation.truncatesins(coeff, max_sins)
            is_truncated = true
        elseif !isnothing(customtruncfunc) && customtruncfunc(mstr, coeff)
            is_truncated = true
        end

        return is_truncated
    end

    msum = mainsum(prop_cache)
    weight_keys = collect(keys(msum.MultiMajoranas))
    empty_weight_keys = [String[] for _ in 1:Threads.maxthreadid()]

    Threads.@threads for i in eachindex(weight_keys)
        tid = Threads.threadid()
        weight_key = weight_keys[i]
        weight_dict = msum.MultiMajoranas[weight_key]

        for ms_int in collect(keys(weight_dict))
            coeff = weight_dict[ms_int]
            if truncfunc(ms_int, coeff)
                delete!(weight_dict, ms_int)
            end
        end

        if isempty(weight_dict)
            push!(empty_weight_keys[tid], weight_key)
        end
    end

    # Remove empty sectors serially to avoid concurrent writes on the top-level Dict.
    for keys_to_delete in empty_weight_keys
        for weight_key in keys_to_delete
            if haskey(msum.MultiMajoranas, weight_key) && isempty(msum.MultiMajoranas[weight_key])
                delete!(msum.MultiMajoranas, weight_key)
            end
        end
    end

    setmainsum!(prop_cache, msum)

    return
end