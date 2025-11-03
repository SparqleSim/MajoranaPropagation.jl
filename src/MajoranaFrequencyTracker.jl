
struct MajoranaFrequencyTracker{CT} <: PathProperties
    coeff::CT
    freq::Int
    nsins::Int
    ncos::Int
end

function MajoranaFrequencyTracker(coeff)
    return MajoranaFrequencyTracker(coeff, 0, 0, 0)
end

function wrapcoefficients(msum::MajoranaSum, ::Type{MProp}) where {MProp<:PathProperties}
    if length(msum) == 0
        throw("The majoranaSum is empty.")
    end

    try
        _, dummy_coeff = first(msum.Majoranas)
        MProp(dummy_coeff)
    catch MethodError
        throw(
            "The constructor `$(MProp)(coeff)` is not defined for the $(MProp) type. " *
            "Either construct a MajoranaSum with wrapped coefficient or define the `$(PProp)(coeff)` constructor.")
    end

    return MajoranaSum(msum.nfermions, Dict(mstr => MProp(coeff) for (mstr, coeff) in msum))
end

function reset_tracker!(msum::MajoranaSum{TT,MajoranaFrequencyTracker{CT}}) where {TT<:Integer,CT}
    for (ms_int, coeff) in msum.Majoranas
        set!(msum, ms_int, MajoranaFrequencyTracker(coeff.coeff, 0, 0, 0))
    end
    return
end

function _applycos(coeff::MajoranaFrequencyTracker{CT}, cos_theta::CT) where {CT}
    return MajoranaFrequencyTracker(coeff.coeff * cos_theta, coeff.freq + 1, coeff.nsins, coeff.ncos + 1)
end
function _applysin(coeff::MajoranaFrequencyTracker{CT}, sin_theta::CT) where {CT}
    return MajoranaFrequencyTracker(coeff.coeff * sin_theta, coeff.freq + 1, coeff.nsins + 1, coeff.ncos)
end
