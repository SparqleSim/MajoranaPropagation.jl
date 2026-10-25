
"""
    VectorMajoranaSum(msum::AbstractMajoranaSum)
    VectorMajoranaSum(mstrs::Union{AbstractArray,Tuple,Base.Generator})

Convert any Majorana sum, or a collection of Majorana terms with coefficients, into a `VectorMajoranaSum`.
"""
VectorMajoranaSum(msum::AbstractMajoranaSum) = VectorMajoranaSum(nsites(msum), is_spinful(msum),
    collect(termtype(msum), terms(msum)), collect(PropagationBase.coefftype(msum), coefficients(msum)))
function VectorMajoranaSum(mstrs::Union{AbstractArray,Tuple,Base.Generator})
    nsites = _checknumberofsites(mstrs)
    is_spinful = _checkspinful(mstrs)

    CType = promote_type(coefftype.(mstrs)...)
    vmsum = VectorMajoranaSum(nsites, is_spinful, [mstr.term for mstr in mstrs], [convert(CType, mstr.coeff) for mstr in mstrs])
    return vmsum
end

"""
    MajoranaSum(msum::AbstractMajoranaSum)

Convert any Majorana sum into a dictionary-based `MajoranaSum`, adding up the coefficients of repeated Majorana strings.
"""
function MajoranaSum(msum::AbstractMajoranaSum)
    # PropagationBase's generic conversion would build a spinless sum
    gathered = MajoranaSum(PropagationBase.coefftype(msum), nsites(msum), is_spinful(msum))
    for (ms, coeff) in msum
        add!(gathered, ms, coeff)
    end
    return gathered
end
