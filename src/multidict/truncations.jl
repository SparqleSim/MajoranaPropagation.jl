function PropagationBase.truncate!(
    prop_cache::MajoranaMultiPropagationCache{MMS};
    max_weight::Real=Inf, min_abs_coeff=1e-10, max_unpaired::Real=Inf,
    max_freq::Real=Inf, max_sins::Real=Inf,
    unpaired_mask=nothing,
    customtruncfunc=nothing,
    kwargs...
) where {TT<:Integer,CT,MMS<:MajoranaSumMulti{TT,CT}}
    if isnothing(unpaired_mask)
        unpaired_mask = create_unpaired_mask(nfermions(mainsum(prop_cache)))
    end

    msum = mainsum(prop_cache)
    weight_keys = collect(keys(msum.MultiMajoranas))
    empty_weight_keys = [Int64[] for _ in 1:Threads.maxthreadid()]
    # Reused thread-local buffers of terms to delete, one per thread.
    to_delete = [TT[] for _ in 1:Threads.maxthreadid()]

    Threads.@threads for i in eachindex(weight_keys)
        tid = Threads.threadid()
        weight_key = weight_keys[i]
        weight_dict = msum.MultiMajoranas[weight_key]
        delete_buf = to_delete[tid]
        empty!(delete_buf)

        for (ms_int, coeff) in weight_dict
            if PauliPropagation.truncatemincoeff(coeff, min_abs_coeff) ||
               truncateunpaired(ms_int, max_unpaired, unpaired_mask) ||
               truncatemajoranaweight(ms_int, max_weight) ||
               PauliPropagation.truncatefrequency(coeff, max_freq) ||
               PauliPropagation.truncatesins(coeff, max_sins) ||
               (!isnothing(customtruncfunc) && customtruncfunc(ms_int, coeff))
                push!(delete_buf, ms_int)
            end
        end

        for ms_int in delete_buf
            if haskey(weight_dict, ms_int)
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