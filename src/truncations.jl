# mask with every odd bit set (bits 1, 3, ..., 2*n_fermions - 1), built by log-doubling;
# pass TT explicitly for a type-stable result (getinttype(n) is only known at runtime)
function create_unpaired_mask(::Type{TT}, n_fermions::Integer) where {TT<:Integer}
    Nbits = 2 * n_fermions
    mask = TT(2)
    s = 2
    while s < Nbits
        mask |= mask << s
        s <<= 1
    end
    full_mask = (TT(1) << Nbits) - TT(1)
    return mask & full_mask
end

create_unpaired_mask(n_fermions::Int) = create_unpaired_mask(getinttype(n_fermions), n_fermions)

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

function PropagationBase.truncate!(
    prop_cache::AbstractMajoranaPropagationCache;
    max_weight::Real=Inf, min_abs_coeff=1e-10, max_unpaired::Real=Inf,
    max_freq::Real=Inf, max_sins::Real=Inf,
    unpaired_mask=nothing,
    customtruncfunc=nothing,
    kwargs...
)
    if isnothing(unpaired_mask) && max_unpaired < Inf
        unpaired_mask = create_unpaired_mask(majoranatype(prop_cache), nfermions(prop_cache))
    end
    function truncfunc(mstr, coeff)
        # slight customization of the truncation function 
        # to truncate majorana weight and single
        is_truncated = false
        if PauliPropagation.truncatemincoeff(coeff, min_abs_coeff)
            is_truncated = true
        elseif max_unpaired < Inf && truncateunpaired(mstr, max_unpaired, unpaired_mask)
            is_truncated = true
        elseif max_weight < Inf && truncatemajoranaweight(mstr, max_weight)
            is_truncated = true
        elseif max_freq < Inf && PauliPropagation.truncatefrequency(coeff, max_freq)
            is_truncated = true
        elseif max_sins < Inf && PauliPropagation.truncatesins(coeff, max_sins)
            is_truncated = true
        elseif !isnothing(customtruncfunc) && customtruncfunc(mstr, coeff)
            is_truncated = true
        end

        return is_truncated
    end
    truncate!(truncfunc, prop_cache; kwargs...)

    return
end

function PropagationBase.truncate!(
    msum::AbstractMajoranaSum;
    max_weight::Real=Inf, min_abs_coeff=1e-10, max_unpaired::Real=Inf,
    max_freq::Real=Inf, max_sins::Real=Inf,
    unpaired_mask=nothing,
    customtruncfunc=nothing,
    kwargs...
)
    if isnothing(unpaired_mask) && max_unpaired < Inf
        unpaired_mask = create_unpaired_mask(majoranatype(msum), nfermions(msum))
    end
    function truncfunc(mstr, coeff)
        # slight customization of the truncation function 
        # to truncate majorana weight and single
        is_truncated = false
        if PauliPropagation.truncatemincoeff(coeff, min_abs_coeff)
            is_truncated = true
        elseif max_unpaired < Inf && truncateunpaired(mstr, max_unpaired, unpaired_mask)
            is_truncated = true
        elseif max_weight < Inf && truncatemajoranaweight(mstr, max_weight)
            is_truncated = true
        elseif max_freq < Inf && PauliPropagation.truncatefrequency(coeff, max_freq)
            is_truncated = true
        elseif max_sins < Inf && PauliPropagation.truncatesins(coeff, max_sins)
            is_truncated = true
        elseif !isnothing(customtruncfunc) && customtruncfunc(mstr, coeff)
            is_truncated = true
        end

        return is_truncated
    end
    msum = truncate!(truncfunc, msum; kwargs...)

    return msum
end