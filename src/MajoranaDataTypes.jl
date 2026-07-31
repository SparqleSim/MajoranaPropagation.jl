using LinearAlgebra
using Bits

"""
    MajoranaString(nfermions::Int, indices::Vector{Int})
    MajoranaString(nfermions::Int, gammas::Int64)

A struct to represent a Majorana string, i.e. a product of Majorana operators, on `nfermions` fermions.
The Majorana operators contained in the string are stored as the set bits of the integer `gammas`, and can alternatively be passed as a vector of `indices`.
A string containing the Majorana indices ``k_1 < k_2 < \\dots < k_w`` represents the Hermitian operator ``i^{\\omega_L} \\, \\gamma_{k_1} \\gamma_{k_2} \\cdots \\gamma_{k_w}`` with ``\\omega_L = w(w-1)/2 \\bmod 2``.
See the `gammas_vector` constructor of `MajoranaSum` for the Majorana index convention.
"""
struct MajoranaString{TT<:Integer}
    nfermions::Int
    gammas::TT
end

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

"""
    nfermions(ms::MajoranaString)

Get the number of fermions that the `MajoranaString` is defined on.
"""
function nfermions(ms::MajoranaString)
    return ms.nfermions
end

### 
abstract type AbstractMajoranaSum <: AbstractTermSum end

majoranas(msum::AbstractMajoranaSum) = terms(msum)

"""
    nfermions(ms::AbstractMajoranaSum)

Get the number of fermions that the `AbstractMajoranaSum` is defined on.
"""
function nfermions(ms::AbstractMajoranaSum)
    if is_spinful(ms)
        return 2 * nsites(ms)
    else
        return nsites(ms)
    end
end

"""
    MajoranaSum{TT<:Integer,CT}

A struct to represent a linear combination of Majorana strings with coefficients of type `CT`, stored as a dictionary mapping the integer representation of each `MajoranaString` to its coefficient.
An entry `(ms, coeff)` represents the operator ``\\mathrm{coeff} \\cdot i^{\\omega_L(ms)} \\, \\gamma_{k_1} \\cdots \\gamma_{k_w}`` with ascending Majorana indices, so every basis term is Hermitian.
See the `gammas_vector` constructor below for the Majorana index convention.
"""
struct MajoranaSum{TT<:Integer,CT} <: AbstractMajoranaSum
    nsites::Int
    is_spinful::Bool
    Majoranas::Dict{TT,CT}
end

# necessary overloads for PropagationBase 
"""
    storage(msum::MajoranaSum)

Get the underlying dictionary of `msum` mapping Majorana string integers to coefficients.
"""
PropagationBase.storage(msum::MajoranaSum) = msum.Majoranas
PropagationBase.nsites(msum::MajoranaSum) = msum.nsites
is_spinful(msum::MajoranaSum) = msum.is_spinful

"""
    MajoranaSum(nfermions::Integer)

Create a MajoranaSum for `nfermions` spinless fermions with `Float64` coefficients.
"""
function MajoranaSum(nfermions::Integer)
    return MajoranaSum(Float64, nfermions)
end

"""
    MajoranaSum(::Type{CT}, n_fermions::Integer) where {CT}

Create a MajoranaSum for `n_fermions` spinless fermions and coefficient type `CT`.
"""
function MajoranaSum(::Type{CT}, n_fermions::Integer) where {CT}
    TT = getinttype(n_fermions)
    is_spinful = false
    return MajoranaSum(n_fermions, is_spinful, Dict{TT,CT}())
end

"""
    MajoranaSum(::Type{CT}, n_sites::Integer, is_spinful::Bool) where {CT}

Create a MajoranaSum with `n_sites` sites, spinful or spinless, and coefficient type `CT`.
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
    MajoranaSum(n_sites::Integer, gammas_vector::Vector{Int}, is_spinful::Bool; coeff=1.)

Create a MajoranaSum with `n_sites` sites and coefficient type `CT` (`Float64` if omitted), containing the single Majorana string given by `gammas_vector` with coefficient `coeff`.
The integers in `gammas_vector` index the Majorana operators that are present:
- spinless fermions: `2 * site - 1` for ``\\gamma`` and `2 * site` for ``\\gamma'``,
- spinful fermions: `4 * site - 3` for ``\\gamma_\\uparrow``, `4 * site - 2` for ``\\gamma'_\\uparrow``, `4 * site - 1` for ``\\gamma_\\downarrow``, `4 * site` for ``\\gamma'_\\downarrow``.
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

"""
    add!(msum::MajoranaSum{TT,CT}, symbol::Symbol, sites, coeff=1.) where {TT<:Integer,CT}
    add!(msum::MajoranaSum{TT,CT}, mstr::MajoranaString{TT}, value::CT) where {TT<:Integer,CT}

Add a term to `msum` in-place: either the observable defined by `symbol` acting on `sites`, scaled by `coeff`, or the single Majorana string `mstr` with coefficient `value`.
See `MajoranaSum(n_sites::Integer, symb::Symbol, sites)` for the supported symbols.
"""
function add!(msum::MajoranaSum{TT,CT}, symbol::Symbol, sites, coeff=1.) where {TT<:Integer,CT}
    add!(msum, coeff * MajoranaSum(nsites(msum), symbol, sites))
    return msum
end

function add!(msum::MajoranaSum{TT,CT}, mstr::MajoranaString{TT}, value::CT) where {TT<:Integer,CT}
    add!(msum, mstr.gammas, value)
end

"""
    set!(ms::MajoranaSum{TT,CT}, ms2::MajoranaString{TT}, value::CT) where {TT<:Integer,CT}

Set the coefficient of the Majorana string `ms2` in `ms` to `value` in-place.
"""
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

"""
    coefftype(::MajoranaSum{TT,CT}) where {TT,CT}

Get the coefficient type `CT` of the `MajoranaSum`.
"""
function coefftype(::MajoranaSum{TT,CT}) where {TT,CT}
    return CT
end

"""
    similar(msum::MajoranaSum)

Create an empty `MajoranaSum` with the same number of sites, spinfulness, and coefficient type as `msum`.
"""
function similar(msum::MajoranaSum)
    new_msum = MajoranaSum(coefftype(msum), nsites(msum), is_spinful(msum))
    sizehint!(new_msum.Majoranas, length(msum.Majoranas))
    return new_msum
end

"""
    get_weight(ms::MajoranaString)
    get_weight(gammas::TT) where {TT<:Integer}

Get the weight of a Majorana string, i.e. the number of Majorana operators it contains.
"""
function get_weight(ms::MajoranaString)
    return get_weight(ms.gammas)
end
function get_weight(gammas::TT) where {TT<:Integer}
    return Bits.weight(gammas)
end

"""
    ==(ms1::MajoranaSum, ms2::MajoranaSum)

Check whether two MajoranaSums are defined on the same system and contain the same Majorana strings with the same coefficients.
"""
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