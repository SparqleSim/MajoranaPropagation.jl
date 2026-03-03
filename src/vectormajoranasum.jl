# vector Majorana sum 

using AcceleratedKernels
const AK = AcceleratedKernels


struct VectorMajoranaSum{TV,CV} <: AbstractMajoranaSum
    nsites::Int
    is_spinful::Bool
    terms::TV
    coeffs::CV

    function VectorMajoranaSum(nsites::Int, is_spinful::Bool, terms::TV, coeffs::CV) where {TV,CV}
        @assert length(terms) == length(coeffs) "Length of terms and coeffs must be the same. Got $(length(terms)) and $(length(coeffs))."
        return new{TV,CV}(nsites, is_spinful, terms, coeffs)
    end
end

# empty initializer for spinless case
VectorMajoranaSum(nsites::Int) = VectorMajoranaSum(Float64, nsites, false)

# empty initializers for both spinless and spinful cases
VectorMajoranaSum(nsites::Int, is_spinful::Bool) = VectorMajoranaSum(Float64, nsites, is_spinful)
VectorMajoranaSum(::Type{CT}, nsites::Int, is_spinful::Bool) where {CT} = VectorMajoranaSum(nsites, is_spinful, getinttype(nsites)[], CT[])

PropagationBase.storage(vmsum::VectorMajoranaSum) = (vmsum.terms, vmsum.coeffs)


"""
    nsites(vmsum::VectorMajoranaSum)

Get the number of qubits that the `VectorMajoranaSum` is defined on.
"""
PropagationBase.nsites(vmsum::VectorMajoranaSum) = vmsum.nsites

""" 
    is_spinful(vmsum::VectorMajoranaSum)
Check if the `VectorMajoranaSum` is defined for spinful fermions.
"""
is_spinful(vmsum::VectorMajoranaSum) = vmsum.is_spinful


Base.similar(vmsum::VectorMajoranaSum) = VectorMajoranaSum(nsites(vmsum), is_spinful(vmsum), Base.similar(vmsum.terms), Base.similar(vmsum.coeffs))

function Base.resize!(vmsum::VectorMajoranaSum, n_new::Int)
    resize!(vmsum.terms, n_new)
    resize!(vmsum.coeffs, n_new)
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
    return vmsum
end