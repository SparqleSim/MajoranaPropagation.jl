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
    truncate!(prop_cache::AbstractMajoranaPropagationCache; max_weight=Inf, min_abs_coeff=1e-10, max_unpaired=Inf, max_freq=Inf, max_sins=Inf, unpaired_mask=nothing, max_terms=nothing, max_terms_rtol=0.05, customtruncfunc=nothing, thread=true, kwargs...)
    truncate!(msum::AbstractMajoranaSum; max_weight=Inf, min_abs_coeff=1e-10, max_unpaired=Inf, max_freq=Inf, max_sins=Inf, unpaired_mask=nothing, max_terms=nothing, max_terms_rtol=0.05, customtruncfunc=nothing, thread=true, kwargs...)

Truncate a Majorana sum (or propagation cache) in-place, removing every Majorana string for which one of the following truncations applies:
- `min_abs_coeff`: the absolute value of the coefficient is smaller than `min_abs_coeff`
- `max_weight`: the weight of the string exceeds `max_weight`
- `max_unpaired`: the number of unpaired Majorana operators (see `compute_unpaired`) exceeds `max_unpaired`; a precomputed mask can be passed as `unpaired_mask`
- `max_freq` / `max_sins`: the frequency / number of sine applications tracked by a `PathProperties` coefficient exceeds the given value
- `customtruncfunc`: a custom function with signature `customtruncfunc(mstr, coeff)::Bool`, returning `true` if the string should be truncated
After these, `max_terms` keeps about this many of the remaining Majorana strings, those with the largest absolute coefficients, within `max_terms_rtol * max_terms` of it (`max_terms_rtol=0` is exact up to ties); see `topkabscoeff`.
Other keyword arguments are ignored.
"""
function PropagationBase.truncate!(mobj::Union{AbstractMajoranaSum,AbstractMajoranaPropagationCache}; thread::Bool=true, kwargs...)
    truncfunc, _, topk = _truncfuncs(mobj; thread, kwargs...)
    truncate!(truncfunc, mobj; thread)
    return _truncatetopk!(mobj, topk; thread)
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
        # a threshold of zero or below truncates nothing; PauliPropagation squares it for complex coefficients
        min_abs_coeff > 0 && PauliPropagation.truncatemincoeff(coeff, min_abs_coeff) && return true
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

# The truncation of `mobj` by the truncation keyword arguments among `kwargs`, as `(truncfunc, termdrop, topk)`:
# the predicate, the part of it that only reads the Majorana string (or `nothing`), and the term budget
# `(max_terms, max_terms_rtol, topk_floor)` of `_truncatetopk!` (or `nothing`), which a gate applies after the predicate,
# once the merged coefficients are known. `topk_floor` is the `Ref` that `propagate!` shares between its gates (or
# `nothing`): it holds half the cutoff of the last cut, below which no string can make the next one, so the predicate
# drops such strings as the gate makes them. Every other keyword argument is ignored, since the gates pass on all of theirs.
function _truncfuncs(mobj; min_abs_coeff::Real=1e-10, max_weight::Real=Inf, max_unpaired::Real=Inf, unpaired_mask=nothing,
    max_freq::Real=Inf, max_sins::Real=Inf, max_terms=nothing, max_terms_rtol::Real=0.05, topk_floor=nothing,
    customtruncfunc=nothing, thread::Bool=true, kwargs...)

    if !isnothing(topk_floor)
        min_abs_coeff = max(min_abs_coeff, topk_floor[])
    end
    truncfunc = buildtruncfunc(mobj; min_abs_coeff, max_weight, max_unpaired, unpaired_mask, max_freq, max_sins, customtruncfunc, thread)
    termdrop = _termtruncfunc(mobj; max_weight, max_unpaired, unpaired_mask)
    topk = isnothing(max_terms) ? nothing : (max_terms, max_terms_rtol, topk_floor)
    return truncfunc, termdrop, topk
end

# Keeps about `k` of the strings of `mobj`, those with the largest absolute coefficients, within `rtol * k` of `k`,
# and leaves the cutoff in `floor` for the next gate (see `_truncfuncs`).
_truncatetopk!(mobj, ::Nothing; thread::Bool=true) = mobj

function _truncatetopk!(mobj, (k, rtol, floor)::Tuple; thread::Bool=true)
    cutoff = topkabscoeff(mobj, k; rtol, thread)
    # zero when `k` is at least the number of strings, so there is nothing to drop
    if cutoff > 0
        truncate!((_, coeff) -> PauliPropagation.truncatemincoeff(coeff, cutoff), mobj; thread)
    end
    if !isnothing(floor)
        # Half the cutoff: a floor at the cutoff itself prunes a little too much (the next gate shrinks the strings just
        # above it, and the cut then cannot bind), half of it keeps the strings between a gate and its cut at about 1.1
        # to 1.5 k and the cut binding. A cut that did not bind halves the floor again instead of lifting it, which would
        # let a gate grow the sum without a floor.
        floor[] = (cutoff > 0 ? cutoff : floor[]) / 2
    end
    return mobj
end


"""
    topkabscoeff(msum::AbstractMajoranaSum, k; rtol=0.05, thread=true)
    topkabscoeff(prop_cache::AbstractMajoranaPropagationCache, k; rtol=0.05, thread=true)

A threshold `t` such that about `k` active coefficients satisfy `abs(coeff) >= t`, within `rtol * k` of `k`.
It is found by radix select on the bits of the magnitudes, never by sorting: a pass counts them into 65536 buckets, and the bucket holding the `k`-th largest is split by further passes on the next 16 bits until keeping or dropping it lands within `rtol * k` of `k`, at most four passes in all.
`rtol=0` is exact up to ties: coefficients of equal magnitude are kept or dropped together, whichever count is closer to `k`.
Returns zero when `k` is at least the number of terms.
"""
function topkabscoeff(thing::Union{AbstractMajoranaSum,AbstractMajoranaPropagationCache}, k::Integer; rtol::Real=0.05, thread::Bool=true)
    k >= 1 || throw(ArgumentError("k must be at least 1, got $k"))
    rtol >= 0 || throw(ArgumentError("rtol must be non-negative, got $rtol"))
    k >= length(thing) && return 0.0

    prefix = UInt64(0)   # the bits above `shift + 16` of the bucket being split
    above = 0            # terms in buckets above the one being split
    for shift in (48, 32, 16, 0)
        counts = _magnitudehistogram(StorageType(thing), thing, prefix, shift, thread)
        for key in (length(counts)-1):-1:0
            with = above + counts[key+1]
            if with < k
                above = with
                continue
            end
            # the whole cutoff bucket is kept or dropped, whichever lands closer to k (ties keep it), once
            # that lands within the tolerance; after the last pass the bucket is a genuine tie
            keep_error = with - k
            drop_error = k - above
            if min(keep_error, drop_error) <= rtol * k || shift == 0
                edge = keep_error <= drop_error ? key : key + 1
                return reinterpret(Float64, ((prefix << 16) + UInt64(edge)) << shift)
            end
            prefix = (prefix << 16) | UInt64(key)
            break
        end
    end
    return 0.0   # unreachable: the bucket being split always holds at least k - above terms
end

# the bits of a non-negative float are monotone in its value
_magnitudebits(coefficient) = reinterpret(UInt64, Float64(abs(PropagationBase.tonumber(coefficient))))

# every storage but a multi sum counts its own coefficients
_magnitudehistogram(::StorageType, thing, prefix::UInt64, shift::Int, thread::Bool) =
    _countmagnitudes(coefficients(thing), prefix, shift, thread)

# a multi sum counts each zone on the zone's own thread, with no tasks inside a zone, and adds the histograms
function _magnitudehistogram(::PropagationBase.MultiSumStorage, thing, prefix::UInt64, shift::Int, thread::Bool)
    count_zone(zone) = _magnitudehistogram(StorageType(zone), zone, prefix, shift, false)
    return sum(PropagationBase._zonevalues(count_zone, Vector{Int}, thing, thread))
end

# Counts the coefficients whose magnitude bits above `shift + 16` equal `prefix` into 2^16 buckets by the
# 16 bits at `shift`. On the first pass `prefix` is 0 and matches everything, since `x >> 64 == 0`.
function _countmagnitudes(coeffs, prefix::UInt64, shift::Int, thread::Bool)
    function count_into!(counts, coefficients)
        @inbounds for coefficient in coefficients
            bits = _magnitudebits(coefficient)
            if bits >> (shift + 16) == prefix
                counts[Int((bits>>shift)&0xFFFF)+1] += 1
            end
        end
        return counts
    end

    # dictionary values are iterated by one task
    if !(coeffs isa AbstractArray)
        return count_into!(zeros(Int, 2^16), coeffs)
    end

    # ponytail: CPU only; a GPU array needs an AK/atomic histogram or a host copy
    PropagationBase._iscpuarray(coeffs) || throw(ArgumentError("topkabscoeff is only supported for CPU arrays, got $(typeof(coeffs))"))
    task_partitioner, n_tasks = PropagationBase._preparetasks(length(coeffs), thread)
    counts = zeros(Int, 2^16, n_tasks)
    count_chunk!(task_id) = count_into!(view(counts, :, task_id), view(coeffs, task_partitioner[task_id]))
    PropagationBase._eachtask(count_chunk!, n_tasks)
    return vec(sum(counts; dims=2))
end
