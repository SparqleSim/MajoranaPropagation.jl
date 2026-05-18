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
        elseif isfinite(max_unpaired) && truncateunpaired(mstr, max_unpaired, unpaired_mask)
            is_truncated = true
        elseif isfinite(max_weight) && truncatemajoranaweight(mstr, max_weight)
            is_truncated = true
        elseif isfinite(max_freq) && PauliPropagation.truncatefrequency(coeff, max_freq)
            is_truncated = true
        elseif isfinite(max_sins) && PauliPropagation.truncatesins(coeff, max_sins)
            is_truncated = true
        elseif !isnothing(customtruncfunc) && customtruncfunc(mstr, coeff)
            is_truncated = true
        end

        return is_truncated
    end

    msum = mainsum(prop_cache)
    @threads for iw in eachindex(msum)
        _truncate_dict!(msum.MultiMajoranas[iw], truncfunc)
    end

    setmainsum!(prop_cache, msum)

    return
end

function PropagationBase.truncate!(
    msum::MajoranaSumMulti;
    max_weight::Real=Inf, min_abs_coeff=1e-10, max_unpaired::Real=Inf,
    max_freq::Real=Inf, max_sins::Real=Inf,
    unpaired_mask=nothing,
    customtruncfunc=nothing,
    kwargs...
)
    if isnothing(unpaired_mask)
        unpaired_mask = create_unpaired_mask(nfermions(msum))
    end
    function truncfunc(mstr, coeff)
        is_truncated = false
        if PauliPropagation.truncatemincoeff(coeff, min_abs_coeff)
            is_truncated = true
        elseif isfinite(max_unpaired) && truncateunpaired(mstr, max_unpaired, unpaired_mask)
            is_truncated = true
        elseif isfinite(max_weight) && truncatemajoranaweight(mstr, max_weight)
            is_truncated = true
        elseif isfinite(max_freq) && PauliPropagation.truncatefrequency(coeff, max_freq)
            is_truncated = true
        elseif isfinite(max_sins) && PauliPropagation.truncatesins(coeff, max_sins)
            is_truncated = true
        elseif !isnothing(customtruncfunc) && customtruncfunc(mstr, coeff)
            is_truncated = true
        end

        return is_truncated
    end

    @threads for iw in eachindex(msum)
        _truncate_dict!(msum.MultiMajoranas[iw], truncfunc)
    end

    return
end

function _truncate_dict!(dict::Dict{TT,CT}, truncfunc::Function) where {TT<:Integer,CT}
    keys_to_delete = TT[]
    for (ms_int, coeff) in dict
        if truncfunc(ms_int, coeff)
            push!(keys_to_delete, ms_int)
        end
    end
    for ms_int in keys_to_delete
        delete!(dict, ms_int)
    end
    return
end