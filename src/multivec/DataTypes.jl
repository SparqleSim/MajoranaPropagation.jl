"""
    MultiVectorMajoranaSum

A weighted-sector Majorana sum where each weight sector is stored as a
`VectorMajoranaSum`. Internally this is a dictionary from Majorana weight
to a `VectorMajoranaSum` restricted to that weight.

This combines:

- the contiguous vector storage of `VectorMajoranaSum`, and
- the parallel weight-sector merging strategy of `MajoranaSumMulti`.
"""
struct MultiVectorMajoranaSum{TT<:Integer,CT} <: AbstractMajoranaSum
    nsites::Int
    is_spinful::Bool
    MultiMajoranas::Dict{Int,VectorMajoranaSum{Vector{TT},Vector{CT}}}
end

"""Construct a `MultiVectorMajoranaSum` from a dense `MajoranaSum`."""
function MultiVectorMajoranaSum(msum::MajoranaSum{TT,CT}) where {TT<:Integer,CT}
    nsites = msum.nsites
    is_spinful = msum.is_spinful

    multimajs = Dict{Int,VectorMajoranaSum{Vector{TT},Vector{CT}}}()

    for (ms_int, coeff) in msum.Majoranas
        weight = get_weight(ms_int)
        vms = get!(multimajs, weight) do
            VectorMajoranaSum(nsites, is_spinful, TT[], CT[])
        end
        push!(vms.terms, ms_int)
        push!(vms.coeffs, coeff)
    end

    return MultiVectorMajoranaSum{TT,CT}(nsites, is_spinful, multimajs)
end

PropagationBase.storage(msum::MultiVectorMajoranaSum) = msum.MultiMajoranas

Base.keys(msum::MultiVectorMajoranaSum) = Base.keys(msum.MultiMajoranas)

function Base.length(msum::MultiVectorMajoranaSum{TT,CT}) where {TT<:Integer,CT}
    total_strings = 0
    for (_w, vms) in msum.MultiMajoranas
        total_strings += length(vms)
    end
    return total_strings
end

function coefftype(::MultiVectorMajoranaSum{TT,CT}) where {TT<:Integer,CT}
    return CT
end

PropagationBase.nsites(msum::MultiVectorMajoranaSum) = msum.nsites

function nfermions(msum::MultiVectorMajoranaSum)
    if msum.is_spinful
        return 2 * msum.nsites
    else
        return msum.nsites
    end
end

"""
    similar(msum::MultiVectorMajoranaSum, W::Int)

Create an empty `MultiVectorMajoranaSum` which has support only on the
sectors with weights `W-2`, `W`, and `W+2`. This mirrors the behaviour
of `similar(::MajoranaSumMulti, W)` and is used by the propagation cache.
"""
#=function similar(msum::MultiVectorMajoranaSum{TT,CT}, W::Int) where {TT<:Integer,CT}
    multimajs = Dict{Int,VectorMajoranaSum{Vector{TT},Vector{CT}}}()
    term_length = max(10, div(length(msum.MultiMajoranas[W]), 2)) # heuristic for preallocating vector sizes in the new sum
    for w in (W - 2, W, W + 2)
        terms = Vector{TT}(undef, term_length)
        coeffs = Vector{CT}(undef, term_length)
        multimajs[w] = VectorMajoranaSum(msum.nsites, msum.is_spinful, terms, coeffs)
    end
    return MultiVectorMajoranaSum{TT,CT}(msum.nsites, msum.is_spinful, multimajs)
end=#

"""
    show_stats(msum::MultiVectorMajoranaSum)

Print statistics on the number of strings per weight sector.
"""
function show_stats(msum::MultiVectorMajoranaSum{TT,CT}) where {TT<:Integer,CT}
    total_strings = length(msum)
    if total_strings == 0
        println("MultiVectorMajoranaSum is empty.")
        return
    end
    sorted_keys = sort(collect(keys(msum.MultiMajoranas)))
    for weight_key in sorted_keys
        nstrings = length(msum.MultiMajoranas[weight_key])
        println("Weight $weight_key: $nstrings strings ($(round(100.0 * nstrings / total_strings))%)")
    end
    println("Total strings: $total_strings")
end