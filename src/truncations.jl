function create_unpaired_mask(n_fermions::Int)
    TT = getinttype(n_fermions)
    mask::TT = 0
    for k = 1:2:(2*n_fermions)
        mask |= TT(1) << k
    end
    return mask
end

function compute_unpaired(res::TT, mask::TT) where {TT<:Integer}
    number_unpaired = res ⊻ (TT(2) * res)
    return Bits.weight(number_unpaired & mask)
end

function truncatemajoranaweight(mstring::MajoranaString, max_weight::Real)
    return get_weight(mstring) > max_weight
end

function truncatemajoranaweight(mstring::TT, max_weight::Real) where {TT<:Integer}
    return get_weight(mstring) > max_weight
end

function truncateunpaired(mstring::TT, max_weight::Real, singles_mask::TT) where {TT<:Integer}
    return compute_unpaired(mstring, singles_mask) > max_weight
end

function create_doublons_filters(Nsites::Int)
    TT = getinttype(2 * Nsites)
    filters::Vector{TT} = []
    for site = 1:Nsites
        filter::TT = 0
        for k = 0:3
            filter |= TT(1) << (4 * (site - 1) + k)
        end
        push!(filters, filter)
    end
    return filters
end

function compute_doublons(res::TT, filters::Vector{TT}) where {TT<:Integer}
    ndoublons = 0
    for filter in filters
        if (res & filter) == filter
            ndoublons += 1
        end
    end
    return ndoublons
end

function checktruncationonall!(
    msum::MajoranaSum{TT,CT}; max_weight::Real=Inf, min_abs_coeff=1e-10, max_unpaired::Real=Inf,
    max_freq::Real=Inf, max_sins::Real=Inf,
    kwargs...
) where {TT<:Integer,CT}
    unpaired_mask = create_unpaired_mask(nfermions(msum))
    for (mstr, coeff) in msum
        mp_checktruncationonone!(
            msum, mstr, coeff, unpaired_mask;
            max_weight=max_weight, min_abs_coeff=min_abs_coeff,
            max_unpaired=max_unpaired,
            max_freq=max_freq, max_sins=max_sins,
            kwargs...
        )
    end
    return
end


function mp_checktruncationonone!(
    msum::MajoranaSum{TT,CT}, mstr::TT, coeff::CT, unpaired_mask::TT;
    max_weight::Real=Inf, min_abs_coeff=1e-10,
    max_unpaired::Real=Inf,
    max_freq::Real=Inf, max_sins::Real=Inf,
    customtruncfunc=nothing,
    kwargs...
) where {TT<:Integer,CT}
    # slight customization of the truncation function
    # to truncate majorana weight and single
    is_truncated = false
    if truncatemajoranaweight(mstr, max_weight)
        is_truncated = true
    elseif truncateunpaired(mstr, max_unpaired, unpaired_mask)
        is_truncated = true
    elseif truncatemincoeff(coeff, min_abs_coeff)
        is_truncated = true
    elseif truncatefrequency(coeff, max_freq)
        is_truncated = true
    elseif truncatesins(coeff, max_sins)
        is_truncated = true
    elseif !isnothing(customtruncfunc) && customtruncfunc(mstr, coeff)
        is_truncated = true
    end
    if is_truncated
        delete!(msum, mstr)
    end
    return
end