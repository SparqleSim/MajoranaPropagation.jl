function create_max_single_filter(N::Int)
    TT = getinttype(N)
    filter::TT = 0
    for k = 1:2:(2*N)
        filter |= TT(1) << k
    end
    return filter
end

function compute_max_single(res::TT, N::Int, filter::TT) where {TT<:Integer}
    number_single = res ⊻ (TT(2) * res)
    return Bits.weight(number_single & filter)
end

function compute_max_single(res::TT, N::Int) where {TT<:Integer}
    number_single = res ⊻ (TT(2) * res)
    filter = create_max_single_filter(N)
    return Bits.weight(number_single & filter)
end

function truncatemajoranaweight(mstring::MajoranaString, max_weight::Real)
    return get_weight(mstring) > max_weight
end

# this overrides the Pauli propagation function if it is called `truncateweight()`
function truncatemajoranaweight(mstring::TT, max_weight::Real) where {TT<:Integer}
    return get_weight(mstring) > max_weight
end

function truncatesingle(mstring::MajoranaString, max_weight::Real)
    return compute_max_single(mstring, mstring.N) > max_weight
end

function truncatesingle(mstring::TT, max_weight::Real, N::Int) where {TT<:Integer}
    return compute_max_single(mstring, N) > max_weight
end

function create_doublons_filters(Nsites::Int)
    TT = getinttype(2 * Nsites)
    filters::Vector{TT} = []
    for site=1:Nsites
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


@inline function PauliPropagation.checktruncationonone!(
    psum::MajoranaSum, pstr, coeff;
    max_weight::Real=Inf, min_abs_coeff=1e-10,
    max_freq::Real=Inf, max_sins::Real=Inf,
    customtruncfunc=nothing,
    kwargs...
)
    # slight customization of the truncation function 
    # to truncate majorana weight and single
    is_truncated = false
    if truncatemajoranaweight(pstr, max_weight)
        is_truncated = true
    #elseif truncatesingle(pstr, maxsingle, psum.nfermions)
    #    is_truncated = true
    elseif truncatemincoeff(coeff, min_abs_coeff)
        is_truncated = true
    elseif truncatefrequency(coeff, max_freq)
        is_truncated = true
    elseif truncatesins(coeff, max_sins)
        is_truncated = true
    elseif !isnothing(customtruncfunc) && customtruncfunc(pstr, coeff)
        is_truncated = true
    end
    if is_truncated
        delete!(psum, pstr)
    end
    return
end