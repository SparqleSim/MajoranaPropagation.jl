# vector Majorana sum 

using AcceleratedKernels
const AK = AcceleratedKernels

const _MIN_ELEMS_PER_TASK = PropagationBase._MIN_ELEMS_PER_TASK

"""
    VectorMajoranaSum(nsites::Int)
    VectorMajoranaSum(nsites::Int, is_spinful::Bool)
    VectorMajoranaSum(::Type{CT}, nsites::Int, is_spinful::Bool) where {CT}

A struct to represent a linear combination of Majorana strings, storing the integer representations and the coefficients in two aligned vectors `terms` and `coeffs`.
"""
mutable struct VectorMajoranaSum{TV,CV} <: AbstractMajoranaSum
    nsites::Int
    is_spinful::Bool
    terms::TV
    coeffs::CV
    # number of leading terms known sorted by integer value and duplicate-free; 0 is always safe
    _terms_sorted::Int

    function VectorMajoranaSum(nsites::Int, is_spinful::Bool, terms::TV, coeffs::CV, _terms_sorted::Int=0) where {TV,CV}
        @assert length(terms) == length(coeffs) "Length of terms and coeffs must be the same. Got $(length(terms)) and $(length(coeffs))."
        @assert 0 <= _terms_sorted <= length(terms) "Sorted prefix length cannot be greater than the number of terms. Got $(length(terms)) and $(_terms_sorted)."
        return new{TV,CV}(nsites, is_spinful, terms, coeffs, _terms_sorted)
    end
end

# empty initializer for spinless case
VectorMajoranaSum(nsites::Int) = VectorMajoranaSum(Float64, nsites, false)

# empty initializers for both spinless and spinful cases
VectorMajoranaSum(nsites::Int, is_spinful::Bool) = VectorMajoranaSum(Float64, nsites, is_spinful)
VectorMajoranaSum(::Type{CT}, nsites::Int, is_spinful::Bool) where {CT} = VectorMajoranaSum(nsites, is_spinful, getinttype(nsites)[], CT[])

"""
    storage(vmsum::VectorMajoranaSum)

Get the tuple `(terms, coeffs)` of the underlying vectors of `vmsum`.
"""
PropagationBase.storage(vmsum::VectorMajoranaSum) = (vmsum.terms, vmsum.coeffs)

PropagationBase.sortedprefix(vmsum::VectorMajoranaSum) = vmsum._terms_sorted
PropagationBase.setsortedprefix!(vmsum::VectorMajoranaSum, n::Int) = (vmsum._terms_sorted = n; vmsum)

majoranatype(vmsum::VectorMajoranaSum{TV,CV}) where {TV,CV} = eltype(TV)


"""
    nsites(vmsum::VectorMajoranaSum)

Get the number of sites that the `VectorMajoranaSum` is defined on.
"""
PropagationBase.nsites(vmsum::VectorMajoranaSum) = vmsum.nsites

"""
    is_spinful(vmsum::VectorMajoranaSum)

Check if the `VectorMajoranaSum` is defined for spinful fermions.
"""
is_spinful(vmsum::VectorMajoranaSum) = vmsum.is_spinful


"""
    similar(vmsum::VectorMajoranaSum)

Create a `VectorMajoranaSum` of the same shape and types as `vmsum`, with uninitialized terms and coefficients.
"""
Base.similar(vmsum::VectorMajoranaSum) = VectorMajoranaSum(nsites(vmsum), is_spinful(vmsum), Base.similar(vmsum.terms), Base.similar(vmsum.coeffs))

"""
    resize!(vmsum::VectorMajoranaSum, n_new::Int)

Resize the terms and coefficients vectors of `vmsum` to length `n_new`.
"""
function Base.resize!(vmsum::VectorMajoranaSum, n_new::Int)
    resize!(vmsum.terms, n_new)
    resize!(vmsum.coeffs, n_new)
    setsortedprefix!(vmsum, min(sortedprefix(vmsum), n_new))  # clamp on shrink, no-op on grow
    return vmsum
end


function Base.show(io::IO, vecmsum::VectorMajoranaSum)
    n_majoranas = length(vecmsum)
    if n_majoranas == 0
        println(io, "Empty VectorMajoranaSum.")
        return
    elseif n_majoranas == 1
        println(io, "VectorMajoranaSum with 1 term:")
    else
        println(io, "VectorMajoranaSum with $(n_majoranas) terms:")
    end

    for i in 1:length(vecmsum)
        if i > 20
            println(io, "  ...")
            break
        end
        println(io, vecmsum.coeffs[i], " * $(reverse(bitstring(vecmsum.terms[i])))")
    end
end


function Base.sort!(vmsum::VectorMajoranaSum; by=nothing, kwargs...)
    # instead of using sortperm, we use sort!() on an index array 
    # this is to be able to sort on any properties of the terms of coeffs 

    indices = collect(1:length(vmsum))

    # default for if "by" is not provided
    byfunc = isnothing(by) ? i -> vmsum.terms[i] : by

    AK.sort!(indices; by=byfunc, kwargs...)
    vmsum.terms .= view(vmsum.terms, indices)
    vmsum.coeffs .= view(vmsum.coeffs, indices)
    setsortedprefix!(vmsum, 0)  # arbitrary order, no dedup: can't assume sorted+unique after this
    return vmsum
end