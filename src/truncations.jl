"""
    create_unpaired_mask(n_fermions::Int)
    create_unpaired_mask(::Type{TT}, n_fermions::Integer) where {TT<:Integer}

Create the bit mask for a system of `n_fermions` fermions to be used in `compute_unpaired`.
Pass the integer type `TT` explicitly for a type-stable result (`getinttype(n_fermions)` is only known at runtime).
"""
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

"""
    compute_unpaired(res::TT, mask::TT) where {TT<:Integer}

Compute the number of unpaired Majorana operators of the Majorana string `res`, i.e. the number of fermionic modes on which exactly one of the two Majorana operators is present.
`mask` is the bit mask created by `create_unpaired_mask`.
"""
function compute_unpaired(res::TT, mask::TT) where {TT<:Integer}
    number_unpaired = res ⊻ (res << 1)
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
    truncate!(prop_cache::AbstractMajoranaPropagationCache; max_weight=Inf, min_abs_coeff=1e-10, max_unpaired=Inf, max_freq=Inf, max_sins=Inf, unpaired_mask=nothing, customtruncfunc=nothing, thread=true, kwargs...)
    truncate!(msum::AbstractMajoranaSum; max_weight=Inf, min_abs_coeff=1e-10, max_unpaired=Inf, max_freq=Inf, max_sins=Inf, unpaired_mask=nothing, customtruncfunc=nothing, thread=true, kwargs...)

Truncate a Majorana sum (or propagation cache) in-place, removing every Majorana string for which one of the following truncations applies:
- `min_abs_coeff`: the absolute value of the coefficient is smaller than `min_abs_coeff`
- `max_weight`: the weight of the string exceeds `max_weight`
- `max_unpaired`: the number of unpaired Majorana operators (see `compute_unpaired`) exceeds `max_unpaired`; a precomputed mask can be passed as `unpaired_mask`
- `max_freq` / `max_sins`: the frequency / number of sine applications tracked by a `PathProperties` coefficient exceeds the given value
- `customtruncfunc`: a custom function with signature `customtruncfunc(mstr, coeff)::Bool`, returning `true` if the string should be truncated
Other keyword arguments are ignored.
"""
function PropagationBase.truncate!(mobj::Union{AbstractMajoranaSum,AbstractMajoranaPropagationCache}; thread::Bool=true, kwargs...)
    truncfunc, _ = _truncfuncs(mobj; thread, kwargs...)
    return truncate!(truncfunc, mobj; thread)
end

"""
    buildtruncfunc(mobj::Union{AbstractMajoranaSum,AbstractMajoranaPropagationCache}; min_abs_coeff=1e-10, max_weight=Inf, max_unpaired=Inf, unpaired_mask=nothing, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, thread=true)

Build the predicate `truncfunc(mstr, coeff)::Bool` that combines the truncations of `truncate!`, returning `true` if the Majorana string `mstr` with coefficient `coeff` should be dropped.
"""
function PauliPropagation.buildtruncfunc(mobj::Union{AbstractMajoranaSum,AbstractMajoranaPropagationCache};
    min_abs_coeff::Real=1e-10, max_weight::Real=Inf, max_unpaired::Real=Inf, unpaired_mask=nothing,
    max_freq::Real=Inf, max_sins::Real=Inf, customtruncfunc=nothing, thread::Bool=true)

    termtrunc = _termtruncfunc(mobj; max_weight, max_unpaired, unpaired_mask)

    function coefftruncfunc(mstr, coeff)
        PauliPropagation.truncatemincoeff(coeff, min_abs_coeff) && return true
        PauliPropagation.truncatefrequency(coeff, max_freq) && return true
        PauliPropagation.truncatesins(coeff, max_sins) && return true
        return !isnothing(customtruncfunc) && customtruncfunc(mstr, coeff)
    end

    # decided once here rather than for every string: without a truncation on the string itself,
    # the predicate never reads the string
    if isnothing(termtrunc)
        return coefftruncfunc
    end

    truncfunc(mstr, coeff) = termtrunc(mstr) || coefftruncfunc(mstr, coeff)
    return truncfunc
end

# The truncation by Majorana weight and by unpaired Majoranas, which only reads the Majorana string, or
# `nothing` if neither is set. A gate may apply it to the strings it creates before it writes them.
function _termtruncfunc(mobj; max_weight::Real=Inf, max_unpaired::Real=Inf, unpaired_mask=nothing)
    if max_weight == Inf && max_unpaired == Inf
        return nothing
    end

    # single assignment to a fresh variable: reassigning the captured kwarg would box it
    mask = (isnothing(unpaired_mask) && max_unpaired < Inf) ?
        create_unpaired_mask(termtype(mobj), nfermions(mobj)) : unpaired_mask

    termtrunc(mstr) = (max_unpaired < Inf && truncateunpaired(mstr, max_unpaired, mask)) ||
                      (max_weight < Inf && truncatemajoranaweight(mstr, max_weight))
    return termtrunc
end

# The truncation of `mobj` by the truncation keyword arguments among `kwargs`, as `(truncfunc, termdrop)`:
# the whole predicate, and the part of it that only reads the Majorana string (or `nothing`).
# Every other keyword argument is ignored, since the gates pass on all of theirs.
function _truncfuncs(mobj; min_abs_coeff::Real=1e-10, max_weight::Real=Inf, max_unpaired::Real=Inf, unpaired_mask=nothing,
    max_freq::Real=Inf, max_sins::Real=Inf, customtruncfunc=nothing, thread::Bool=true, kwargs...)

    truncfunc = buildtruncfunc(mobj; min_abs_coeff, max_weight, max_unpaired, unpaired_mask, max_freq, max_sins, customtruncfunc, thread)
    termdrop = _termtruncfunc(mobj; max_weight, max_unpaired, unpaired_mask)
    return truncfunc, termdrop
end
