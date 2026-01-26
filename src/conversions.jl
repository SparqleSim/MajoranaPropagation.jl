
VectorMajoranaSum(msum::MajoranaSum) = VectorMajoranaSum(msum.nsites, msum.is_spinful, collect(majoranas(msum)), collect(coefficients(msum)))
function VectorMajoranaSum(mstrs::Union{AbstractArray,Tuple,Base.Generator})
    nsites = _checknumberofsites(mstrs)
    is_spinful = _checkspinful(mstrs)

    CType = promote_type(coefftype.(mstrs)...)
    vmsum = VectorMajoranaSum(nsites, is_spinful, [mstr.term for mstr in mstrs], [convert(CType, mstr.coeff) for mstr in mstrs])
    return vmsum
end


#=# Conversion of MajoranaString and MajoranaSum to different coefficient types
function Base.convert(::Type{MajoranaString{TT1,CT1}}, mstr::MajoranaString{TT2,CT2}) where {TT1,TT2,CT1,CT2}
    if TT1 != TT2
        throw(ArgumentError("Cannot change term type from $TT2 to $TT1"))
    end
    return MajoranaString(mstr.nsites, mstr.is_spinful, convert(TT1, mstr.term), convert(CT1, mstr.coeff))
end

function Base.convert(::Type{MajoranaSum{TT1,CT1}}, msum::MajoranaSum{TT2,CT2}) where {TT1,TT2,CT1,CT2}
    if TT1 != TT2
        throw(ArgumentError("Cannot change term type from $TT2 to $TT1"))
    end
    return MajoranaSum(msum.nsites, msum.is_spinful, convert(Dict{TT1,CT1}, msum.terms))
end

function convertcoefftype(::Type{CT1}, mstr::MajoranaString{TT,CT2}) where {TT,CT1,CT2}
    return MajoranaString(mstr.nsites, mstr.is_spinful, mstr.term, convert(CT1, mstr.coeff))
end

function convertcoefftype(::Type{CT1}, msum::MajoranaSum{TT,CT2}) where {TT,CT1,CT2}
    return MajoranaSum(msum.nsites, msum.is_spinful, convert(Dict{TT,CT1}, msum.terms))
end


# Checks whether the number of sites is the same between datatypes.
function _checknumberofsites(nsites::Integer, mobj)
    obj_nsites = nsites(mobj)
    if nsites != obj_nsites
        throw(
            ArgumentError(
                "Number of sites ($(nsites)) must equal number of sites ($(obj_nsites)) in $(typeof(mobj))"
            )
        )
    end
    return nsites
end

function _checknumberofsites(mobj1, mobj2)
    nsites1 = nsites(mobj1)
    nsites2 = nsites(mobj2)
    if nsites1 != nsites2
        throw(
            ArgumentError(
                "Number of sites ($(nsites1)) in $(typeof(mobj1)) must equal number of sites ($(nsites2)) in $(typeof(mobj2))"
            )
        )
    end
    return nsites1
end

function _checknumberofsites(mobjects::Union{AbstractArray,Tuple,Base.Generator})
    if !allequal(nsites(mobj) for mobj in mobjects)
        throw(
            ArgumentError(
                "Number of sites in passed collection of type $(typeof(mobjects)) is not consistent."
            )
        )
    end
    return nsites(first(mobjects))
end

# Checks whether spinful property is consistent.
function _checkspinful(is_spinful::Bool, mobj)
    obj_spinful = is_spinful(mobj)
    if is_spinful != obj_spinful
        throw(
            ArgumentError(
                "Spinful property ($(is_spinful)) must equal spinful property ($(obj_spinful)) in $(typeof(mobj))"
            )
        )
    end
    return is_spinful
end

function _checkspinful(mobj1, mobj2)
    spinful1 = is_spinful(mobj1)
    spinful2 = is_spinful(mobj2)
    if spinful1 != spinful2
        throw(
            ArgumentError(
                "Spinful property ($(spinful1)) in $(typeof(mobj1)) must equal spinful property ($(spinful2)) in $(typeof(mobj2))"
            )
        )
    end
    return spinful1
end

function _checkspinful(mobjects::Union{AbstractArray,Tuple,Base.Generator})
    if !allequal(is_spinful(mobj) for mobj in mobjects)
        throw(
            ArgumentError(
                "Spinful property in passed collection of type $(typeof(mobjects)) is not consistent."
            )
        )
    end
    return is_spinful(first(mobjects))
end

# throw error for mismatched term types
function _checktermtype(mobj1, mobj2)
    if majoranatype(mobj1) != majoranatype(mobj2)
        throw(ArgumentError("Majorana types do not match. Got $(majoranatype(mobj1)) and $(majoranatype(mobj2))."))
    end
end=#