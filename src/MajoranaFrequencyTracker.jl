
"""
    MajoranaFrequencyTracker(coeff)

`PathProperties` coefficient wrapper tracking the number of splits (`freq`), sine (`nsins`)
and cosine (`ncos`) applications along a Majorana propagation path.
"""
struct MajoranaFrequencyTracker{CT} <: PathProperties
    coeff::CT
    freq::Int
    nsins::Int
    ncos::Int
end

MajoranaFrequencyTracker(coeff::Number) = MajoranaFrequencyTracker(float(coeff), 0, 0, 0)
MajoranaFrequencyTracker{CT}(coeff::Number) where {CT} = MajoranaFrequencyTracker{CT}(convert(CT, coeff), 0, 0, 0)

PropagationBase.numcoefftype(::Type{MajoranaFrequencyTracker{CT}}) where {CT} = CT

"""
    wrapcoefficients(msum, ::Type{MProp}) where {MProp<:PathProperties}

Wrap the coefficients of a `MajoranaSum` or `VectorMajoranaSum` into the `PathProperties`
type `MProp` via its one-argument constructor `MProp(coeff)`.

For `VectorMajoranaSum` the sorted prefix is preserved. Note that path counters are only
exact on the vector backend while merges take the sorted-tail path: the generic fallback
merge (unsorted prefix, gate weight > 4) resets the counters of merged terms.
"""
function wrapcoefficients(msum::MajoranaSum, ::Type{MProp}) where {MProp<:PathProperties}
    return MajoranaSum(msum.nsites, msum.is_spinful, Dict(mstr => MProp(coeff) for (mstr, coeff) in msum.Majoranas))
end

function wrapcoefficients(vmsum::VectorMajoranaSum, ::Type{MProp}) where {MProp<:PathProperties}
    return VectorMajoranaSum(vmsum.nsites, vmsum.is_spinful, copy(vmsum.terms), map(MProp, vmsum.coeffs), vmsum._terms_sorted)
end

"""
    unwrapcoefficients(msum)

Return a copy of the sum with the `PathProperties` wrappers replaced by their `coeff` field.
"""
function unwrapcoefficients(msum::MajoranaSum{TT,<:PathProperties}) where {TT}
    return MajoranaSum(msum.nsites, msum.is_spinful, Dict(mstr => coeff.coeff for (mstr, coeff) in msum.Majoranas))
end

function unwrapcoefficients(vmsum::VectorMajoranaSum{TV,<:AbstractVector{<:PathProperties}}) where {TV}
    return VectorMajoranaSum(vmsum.nsites, vmsum.is_spinful, copy(vmsum.terms), [coeff.coeff for coeff in vmsum.coeffs], vmsum._terms_sorted)
end

"""
    reset_tracker!(msum)

Reset the `freq`, `nsins` and `ncos` counters of all coefficients to 0, keeping the coefficients.
"""
function reset_tracker!(msum::MajoranaSum{TT,MajoranaFrequencyTracker{CT}}) where {TT<:Integer,CT}
    for (ms_int, coeff) in msum.Majoranas
        set!(msum, ms_int, MajoranaFrequencyTracker(coeff.coeff, 0, 0, 0))
    end
    return
end

function reset_tracker!(vmsum::VectorMajoranaSum{TV,<:AbstractVector{<:MajoranaFrequencyTracker}}) where {TV}
    map!(coeff -> MajoranaFrequencyTracker(coeff.coeff, 0, 0, 0), vmsum.coeffs, vmsum.coeffs)
    return
end

function _applycos(coeff::MajoranaFrequencyTracker, cos_theta::Number)
    return MajoranaFrequencyTracker(coeff.coeff * cos_theta, coeff.freq + 1, coeff.nsins, coeff.ncos + 1)
end
function _applysin(coeff::MajoranaFrequencyTracker, sin_theta::Number)
    return MajoranaFrequencyTracker(coeff.coeff * sin_theta, coeff.freq + 1, coeff.nsins + 1, coeff.ncos)
end
