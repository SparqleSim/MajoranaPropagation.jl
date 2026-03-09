using LinearAlgebra
using Bits

struct MajoranaString{TT<:Integer}
    nfermions::Int
    gammas::TT
end

# TODO: documentation
function MajoranaString(nfermions::Int, indices::Vector{Int})
    TT = getinttype(nfermions)
    gammas = _bitonesat(TT, indices)
    return MajoranaString(nfermions, gammas)
end

function MajoranaString(nfermions::Int, gammas::Int64)
    # Int64 is probably unwanted, lets make it the correct type
    TT = getinttype(nfermions)
    return MajoranaString(nfermions, convert(TT, gammas))
end

function nfermions(ms::MajoranaString)
    return ms.nfermions
end

### 
abstract type AbstractMajoranaSum <: AbstractTermSum end

majoranas(msum::AbstractMajoranaSum) = terms(msum)

""" 
    nfermions(vmsum::AbstractMajoranaSum)

    Get the number of fermions that the `AbstractMajoranaSum` is defined on.
"""
function nfermions(ms::AbstractMajoranaSum)
    if is_spinful(ms)
        return 2 * nsites(ms)
    else
        return nsites(ms)
    end
end

struct MajoranaSum{TT<:Integer,CT} <: AbstractMajoranaSum
    nsites::Int
    is_spinful::Bool
    Majoranas::Dict{TT,CT}
end

# necessary overloads for PropagationBase 
PropagationBase.storage(msum::MajoranaSum) = msum.Majoranas
PropagationBase.nsites(msum::MajoranaSum) = msum.nsites
is_spinful(msum::MajoranaSum) = msum.is_spinful

""" 
    MajoranaSum(n_fermions::Integer)
Create a MajoranaSum for `nfermions` spinless fermions with Float64 coefficients.
"""
function MajoranaSum(nfermions::Integer)
    return MajoranaSum(Float64, nfermions)
end

""" 
    MajoranaSum(::Type{CT}, n_fermions::Integer) where {CT}
Create a MajoranaSum for `nfermions` spinless fermions and coefficient type `CT`.
"""
function MajoranaSum(::Type{CT}, n_fermions::Integer) where {CT}
    TT = getinttype(n_fermions)
    is_spinful = false
    return MajoranaSum(n_fermions, is_spinful, Dict{TT,CT}())
end

""" 
    MajoranaSum(::Type{CT}, n_sites::Integer, is_spinful::Bool) where {CT}
Create a MajoranaSum for with `n_sites` that can be both spinful or spinless (depending on `is_spinful::Bool`) and coefficient type `CT`.
"""
function MajoranaSum(::Type{CT}, n_sites::Integer, is_spinful::Bool) where {CT}
    if is_spinful
        TT = getinttype(2 * n_sites)
    else
        TT = getinttype(n_sites)
    end
    return MajoranaSum(n_sites, is_spinful, Dict{TT,CT}())
end

"""
    MajoranaSum(::Type{CT}, n_sites::Integer, gammas_vector::Vector{Int}, is_spinful::Bool; coeff=1.) where {CT}
Create a MajoranaSum for with `n_sites` and coefficient type `CT`
with a specific initial configuration of Majorana operators given by a list of integers indicating which Majorana operators are present
which is indexed as 
    [spinless fermions]:
        index = 2 * site -1 for gamma
        index = 2 * site for gamma prime
and 
    [spinful fermions]:
        index = 4 * site - 3 for gamma up
        index = 4 * site - 2  for gamma prime up
        index = 4 * site - 1 for gamma down
        index = 4 * site for gamma prime down
"""
function MajoranaSum(::Type{CT}, n_sites::Integer, gammas_vector::Vector{Int}, is_spinful::Bool; coeff=1.) where {CT}
    coeff = CT(coeff)
    n_fermions = is_spinful ? 2 * n_sites : n_sites
    mstring = MajoranaString(n_fermions, gammas_vector)
    return MajoranaSum(n_sites, is_spinful, Dict(mstring.gammas => coeff))
end

function MajoranaSum(n_sites::Integer, gammas_vector::Vector{Int}, is_spinful::Bool; coeff=1.)
    return MajoranaSum(Float64, n_sites, gammas_vector, is_spinful; coeff=coeff)
end


import PauliPropagation.PropagationBase: add!, set!, delete!, empty!

function add!(ms::MajoranaSum{TT,CT}, symbol::Symbol, sites, coeff=1.) where {TT<:Integer,CT}
    add!(ms, coeff * MajoranaSum(nsites(ms), symbol, sites))
    return ms
end

function add!(ms::MajoranaSum{TT,CT}, ms2::MajoranaString{TT}, value::CT) where {TT<:Integer,CT}
    add!(ms, ms2.gammas, value)
end

function set!(ms::MajoranaSum{TT,CT}, ms2::MajoranaString{TT}, value::CT) where {TT<:Integer,CT}
    set!(ms, ms2.gammas, value)
    return
end

function Base.pop!(ms::MajoranaSum{TT,CT}, ms2_gammas::TT) where {TT<:Integer,CT}
    return pop!(ms.Majoranas, ms2_gammas, 0.)
end


function Base.mergewith!(merge, msum1::MajoranaSum, msum2::MajoranaSum)
    mergewith!(merge, msum1.Majoranas, msum2.Majoranas)
    return msum1
end

function Base.show(io::IO, ms::MajoranaString)
    print(io, "$(reverse(bitstring(ms.gammas)))")
end

function Base.show(io::IO, ms::MajoranaSum)
    max_display = 20
    print(io, "MajoranaSum with $(length(ms)) terms:")
    for (i, (mstring, coeff)) in enumerate(ms.Majoranas)
        if i <= max_display
            print(io, "\n")
            print(io, "    $(coeff) * $(reverse(bitstring(mstring)))")
        else
            print(io, "\n    ...")
            break
        end
    end
end


function majoranatype(::MajoranaSum{TT,CT}) where {TT,CT}
    return TT
end

function coefftype(::MajoranaSum{TT,CT}) where {TT,CT}
    return CT
end

function similar(msum::MajoranaSum)
    new_msum = MajoranaSum(coefftype(msum), nsites(msum), is_spinful(msum))
    sizehint!(new_msum.Majoranas, length(msum.Majoranas))
    return new_msum
end

function get_weight(ms::MajoranaString)
    return get_weight(ms.gammas)
end
function get_weight(gammas::TT) where {TT<:Integer}
    return Bits.weight(gammas)
end

function Base.:(==)(ms1::MajoranaSum, ms2::MajoranaSum)
    if nsites(ms1) != nsites(ms2)
        return false
    end
    if is_spinful(ms1) != is_spinful(ms2)
        return false
    end
    return ms1.Majoranas == ms2.Majoranas
end

function pop_id!(msum::MajoranaSum)
    if haskey(msum.Majoranas, 0)
        delete!(msum.Majoranas, 0)
    end
    return
end

# a function to get bits=1 at specified positions
# indices here is some sort of iterable
function _bitonesat(::Type{TT}, indices) where {TT<:Integer}
    mask = zero(TT)
    for pos in indices
        mask |= TT(1) << (pos - 1)
    end
    return mask
end

function _bitonesat(::Type{TT}, index::Integer) where {TT<:Integer}
    return TT(1) << (index - 1)
end

function _checknfermions(msum1::MajoranaSum, msum2::MajoranaSum)
    if nfermions(msum1) != nfermions(msum2)
        throw(ArgumentError("MajoranaSums must have the same nfermions, but have $(nfermions(msum1)) and $(nfermions(msum2))"))
    end

end

function _checknfermions(msum::MajoranaSum, ms::MajoranaString)
    if nfermions(msum) != nfermions(ms)
        throw(ArgumentError("MajoranaSum and MajoranaString must have the same nfermions, but have $(nfermions(msum)) and $(nfermions(ms))"))
    end
end

function _checknfermions(ms1::MajoranaString, ms2::MajoranaString)
    if nfermions(ms1) != nfermions(ms2)
        throw(ArgumentError("Majorana strings must have the same length, but have lengths $(nfermions(ms1)) and $(nfermions(ms2))"))
    end
end


include("vectormajoranasum.jl")
include("conversions.jl")