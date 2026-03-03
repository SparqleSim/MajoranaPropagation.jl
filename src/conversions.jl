
VectorMajoranaSum(msum::MajoranaSum) = VectorMajoranaSum(msum.nsites, msum.is_spinful, collect(majoranas(msum)), collect(coefficients(msum)))
function VectorMajoranaSum(mstrs::Union{AbstractArray,Tuple,Base.Generator})
    nsites = _checknumberofsites(mstrs)
    is_spinful = _checkspinful(mstrs)

    CType = promote_type(coefftype.(mstrs)...)
    vmsum = VectorMajoranaSum(nsites, is_spinful, [mstr.term for mstr in mstrs], [convert(CType, mstr.coeff) for mstr in mstrs])
    return vmsum
end
