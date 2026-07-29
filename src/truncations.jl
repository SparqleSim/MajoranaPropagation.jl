"""
    create_unpaired_mask(n_fermions::Int)

Create the bit mask for a system of `n_fermions` fermions to be used in `compute_unpaired`.
"""
function create_unpaired_mask(n_fermions::Int)
    TT = getinttype(n_fermions)
    mask::TT = 0
    for k = 1:2:(2*n_fermions)
        mask |= TT(1) << k
    end
    return mask
end

"""
    compute_unpaired(res::TT, mask::TT) where {TT<:Integer}

Compute the number of unpaired Majorana operators of the Majorana string `res`, i.e. the number of fermionic modes on which exactly one of the two Majorana operators is present.
`mask` is the bit mask created by `create_unpaired_mask`.
"""
function compute_unpaired(res::TT, mask::TT) where {TT<:Integer}
    number_unpaired = res ⊻ (TT(2) * res)
    return Bits.weight(number_unpaired & mask)
end

"""
    truncatemajoranaweight(mstring::MajoranaString, max_weight::Real)
    truncatemajoranaweight(mstring::TT, max_weight::Real) where {TT<:Integer}

Return `true` if the weight of the Majorana string `mstring` exceeds `max_weight`.
"""
function truncatemajoranaweight(mstring::MajoranaString, max_weight::Real)
    return get_weight(mstring) > max_weight
end

function truncatemajoranaweight(mstring::TT, max_weight::Real) where {TT<:Integer}
    return get_weight(mstring) > max_weight
end

function truncateunpaired(mstring::TT, max_weight::Real, singles_mask::TT) where {TT<:Integer}
    return compute_unpaired(mstring, singles_mask) > max_weight
end

"""
    create_doublons_filters(Nsites::Int)

Create one bit filter per site of a spinful system, each covering the four Majorana operators of that site, to be used in `compute_doublons`.
"""
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

"""
    compute_doublons(res::TT, filters::Vector{TT}) where {TT<:Integer}

Compute the number of sites of a spinful system on which the Majorana string `res` contains all four Majorana operators, using the `filters` created by `create_doublons_filters`.
"""
function compute_doublons(res::TT, filters::Vector{TT}) where {TT<:Integer}
    ndoublons = 0
    for filter in filters
        if (res & filter) == filter
            ndoublons += 1
        end
    end
    return ndoublons
end

"""
    truncate!(prop_cache::AbstractMajoranaPropagationCache; max_weight=Inf, min_abs_coeff=1e-10, max_unpaired=Inf, max_freq=Inf, max_sins=Inf, unpaired_mask=nothing, customtruncfunc=nothing, kwargs...)
    truncate!(msum::AbstractMajoranaSum; max_weight=Inf, min_abs_coeff=1e-10, max_unpaired=Inf, max_freq=Inf, max_sins=Inf, unpaired_mask=nothing, customtruncfunc=nothing, kwargs...)

Truncate a Majorana sum (or propagation cache) in-place, removing every Majorana string for which one of the following truncations applies:
- `min_abs_coeff`: the absolute value of the coefficient is smaller than `min_abs_coeff`
- `max_weight`: the weight of the string exceeds `max_weight`
- `max_unpaired`: the number of unpaired Majorana operators (see `compute_unpaired`) exceeds `max_unpaired`; a precomputed mask can be passed as `unpaired_mask`
- `max_freq` / `max_sins`: the frequency / number of sine applications tracked by a `PathProperties` coefficient exceeds the given value
- `customtruncfunc`: a custom function with signature `customtruncfunc(mstr, coeff)::Bool`, returning `true` if the string should be truncated
"""
function PropagationBase.truncate!(
    prop_cache::AbstractMajoranaPropagationCache;
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
        # slight customization of the truncation function 
        # to truncate majorana weight and single
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
    if isnothing(unpaired_mask)
        unpaired_mask = create_unpaired_mask(nfermions(msum))
    end
    function truncfunc(mstr, coeff)
        # slight customization of the truncation function 
        # to truncate majorana weight and single
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
    msum = truncate!(truncfunc, msum; kwargs...)

    return msum
end