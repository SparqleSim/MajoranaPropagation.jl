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

struct MajoranaSum{TT<:Integer,CT}
    nsites::Int
    is_spinful::Bool
    Majoranas::Dict{TT,CT}
end

""" 
    MajoranaSum(n_fermions::Integer)
Create a MajoranaSum for `nfermions` spinless fermions and coefficient type `CT`.
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

function add!(msum::MajoranaSum{TT,CT}, symbol::Symbol, sites) where {TT<:Integer,CT}
    add!(msum, MajoranaSum(msum.nsites, symbol, sites))
    return msum
end

function add!(msum::MajoranaSum{TT,CT}, ms2::MajoranaSum{TT,CT}) where {TT<:Integer,CT}
    mergewith!(+, msum.Majoranas, ms2.Majoranas)
    return msum
end


function add!(msum::MajoranaSum{TT,CT}, ms::MajoranaString{TT}, value::CT) where {TT<:Integer,CT}
    add!(msum, ms.gammas, value)
end

function add!(msum::MajoranaSum, ms_gammas::TT, value::CT) where {TT<:Integer,CT}
    if haskey(msum.Majoranas, ms_gammas)
        msum.Majoranas[ms_gammas] += value
    else
        msum.Majoranas[ms_gammas] = value
    end
end

function set!(msum::MajoranaSum{TT,CT}, ms::MajoranaString, value::CT) where {TT<:Integer,CT}
    set!(msum, ms.gammas, value)
    return
end

function set!(msum::MajoranaSum{TT,CT}, ms::TT, value::CT) where {TT<:Integer,CT}
    msum.Majoranas[ms] = value
    return
end

function majoranas(msum::MajoranaSum)
    return keys(msum.Majoranas)
end

function coefficients(msum::MajoranaSum)
    return values(msum.Majoranas)
end

function nfermions(msum::MajoranaSum)
    if msum.is_spinful
        return 2 * msum.nsites
    else
        return msum.nsites
    end
end

function nfermions(ms::MajoranaString)
    return ms.nfermions
end

function Base.delete!(msum::MajoranaSum{TT,CT}, ms::MajoranaString{TT}) where {TT<:Integer,CT}
    delete!(msum.Majoranas, ms.gammas)
end
function Base.delete!(msum::MajoranaSum{TT,CT}, ms_gammas::TT) where {TT<:Integer,CT}
    delete!(msum.Majoranas, ms_gammas)
end

function Base.pop!(msum::MajoranaSum{TT,CT}, ms_gammas::TT) where {TT<:Integer,CT}
    return pop!(msum.Majoranas, ms_gammas, 0.)
end

function Base.length(msum::MajoranaSum)
    return length(msum.Majoranas)
end

function Base.mergewith!(merge, msum1::MajoranaSum, msum2::MajoranaSum)
    mergewith!(merge, msum1.Majoranas, msum2.Majoranas)
    return msum1
end

function Base.empty!(msum::MajoranaSum)
    empty!(msum.Majoranas)
    return msum
end

function Base.show(io::IO, ms::MajoranaString)
    print(io, "$(reverse(string(ms.gammas; base=2, pad=2 * ms.nfermions)))")
end

function Base.show(io::IO, msum::MajoranaSum)
    max_display = 8
    print(io, "MajoranaSum with $(length(msum)) term$(length(msum) == 1 ? "" : "s"):(")
    for (i, (mstring, coeff)) in enumerate(msum.Majoranas)
        if i <= max_display
            print(io, "\n")
            print(io, "    $(coeff) * $(reverse(string(mstring; base=2, pad=2 * nfermions(msum))))")
        else
            print(io, "\n    ...")
            break
        end
    end
    print(io, ")")
end


function majoranatype(::MajoranaSum{TT,CT}) where {TT,CT}
    return TT
end

function coefftype(::MajoranaSum{TT,CT}) where {TT,CT}
    return CT
end

function similar(msum::MajoranaSum)
    new_msum = MajoranaSum(coefftype(msum), msum.nsites, msum.is_spinful)
    sizehint!(new_msum.Majoranas, length(msum.Majoranas))
    return new_msum
end

Base.iterate(msum::MajoranaSum, state=1) = iterate(msum.Majoranas, state)

function get_weight(ms::MajoranaString)
    return get_weight(ms.gammas)
end
function get_weight(gammas::TT) where {TT<:Integer}
    return Bits.weight(gammas)
end

function compute_parity_bits_and_shift(u::TT, Nbits::Int) where {TT<:Integer}

    # If Nbits=1 there is no parity
    if Nbits <= 1
        return TT(0)
    end

    # TODO: these masks can be precomputed for efficiency

    # mask for all active bits
    full_mask = (TT(1) << Nbits) - TT(1)

    # mask for Nbits - 1 bits.
    mask = (full_mask >> 1)

    # crop last bit
    p = u & mask

    # this is a parallel prefix xor operation
    # runs in log2(Nbits) steps
    s = 1
    while s < Nbits
        p ⊻= (p << s)
        s <<= 1
    end

    # shift necessary for consistency with site convention
    p = p << 1

    # mask all bits
    return p & full_mask
end

function omega_L_mult(ms1::MajoranaString, ms2::MajoranaString)
    return omega_L_mult(ms1.gammas, ms2.gammas, 2 * ms1.nfermions)
end

function omega_L_mult(ms1::TT, ms2::TT, Nbits) where {TT<:Integer}
    return mod(Bits.weight(ms1 & compute_parity_bits_and_shift(ms2, Nbits)), 2)
end

function omega_L_mult(ms::TT) where {TT<:Integer}
    wms = get_weight(ms)
    return mod((wms^2 - wms) / 2, 2)
end

function omega_L_mult(ms::MajoranaString)
    return omega_L_mult(ms.gammas)
end

function omega_mult(ms1::MajoranaString, ms2::MajoranaString)
    return omega_mult(ms1.gammas, ms2.gammas)
end

function omega_mult(gammas1::TT, gammas2::TT) where {TT<:Integer}
    w1 = get_weight(gammas1)
    w2 = get_weight(gammas2)
    return mod(w1 * w2 - get_weight(gammas1 & gammas2), 2)
end

function omega_mult(ms::MajoranaString)
    return omega_L_mult(ms, ms)
end

function Base.:(==)(ms1::MajoranaSum, ms2::MajoranaSum)
    if ms1.nsites != ms2.nsites
        return false
    end
    if ms1.is_spinful != ms2.is_spinful
        return false
    end
    return ms1.Majoranas == ms2.Majoranas
end

function mstring_additon(ms1::TT, ms2::TT) where {TT<:Integer}
    return ms1 ⊻ ms2
end
function Base.:(+)(ms1::MajoranaString, ms2::MajoranaString)
    _checknfermions(ms1, ms2)
    return MajoranaString(ms1.nfermions, mstring_additon(ms1.gammas, ms2.gammas))
end


function Base.:(+)(msum1::MajoranaSum, msum2::MajoranaSum)
    _checknfermions(msum1, msum2)
    msum1 = deepcopy(msum1)
    add!(msum1, msum2)
    return msum1
end

function Base.:(*)(msum1::MajoranaSum, msum2::MajoranaSum)
    _checknfermions(msum1, msum2)
    res = MajoranaSum(coefftype(msum1), msum1.nsites, msum1.is_spinful)
    for (ms1, coeff1) in msum1.Majoranas
        for (ms2, coeff2) in msum2.Majoranas
            prefactor, ms3 = ms_mult(ms1, ms2, nfermions(msum1))
            @assert imag(prefactor) ≈ 0
            prefactor = real(prefactor)
            add!(res, ms3, prefactor * tonumber(coeff1) * tonumber(coeff2))
        end
    end
    return res
end

function Base.:(*)(coeff::CT, msum::MajoranaSum{TT,CT}) where {TT<:Integer,CT}
    res = similar(msum)
    for (ms1, coeff1) in msum.Majoranas
        set!(res, ms1, coeff * coeff1)
    end
    return res
end

function fprefactor(g1::TT, g2::TT) where {TT<:Integer}
    return omega_L_mult(g1) * omega_L_mult(g2) + omega_mult(g1, g2) * (omega_L_mult(g1) + omega_L_mult(g2) + 1)
end

function fprefactor(ms1::MajoranaString, ms2::MajoranaString)
    return fprefactor(ms1.gammas, ms2.gammas)
end

function ms_mult(ms1::MajoranaString, ms2::MajoranaString)
    if ms1.nfermions != ms2.nfermions
        throw(ArgumentError("Majorana strings must have the same length, but have lengths $(ms1.nfermions) and $(ms2.nfermions)"))
    end
    prefactor, result = ms_mult(ms1.gammas, ms2.gammas, 2 * ms1.nfermions)
    return prefactor, MajoranaString(ms1.nfermions, result)
end

function ms_mult(ms1::TT, ms2::TT, n_fermions::Integer) where {TT<:Integer}
    result = mstring_additon(ms1, ms2) # result = ms1 + ms2
    prefactor = (-1)^(omega_L_mult(ms1, ms2, 2 * n_fermions) + fprefactor(ms1, ms2))
    if mod(omega_mult(ms1, ms2), 2) == 1
        return 1im * prefactor, result
    end
    return prefactor, result
end

function commutes(ms1::MajoranaString, ms2::MajoranaString)
    return commutes(ms1.gammas, ms2.gammas)
end

function commutes(gammas1::Integer, gammas2::Integer)
    return mod(omega_mult(gammas1, gammas2), 2) == 0
end


function norm(msum::MajoranaSum, L=2)
    if length(msum) == 0
        return 0.0
    end
    return LinearAlgebra.norm((coeff for coeff in coefficients(msum)), L)
end

function commutator(msum1::MajoranaSum, msum2::MajoranaSum)
    res = MajoranaSum(msum1.nfermions, typeof(1.1im))
    for (ms1, coeff1) in msum1.Majoranas
        for (ms2, coeff2) in msum2.Majoranas
            if commutes(ms1, ms2)
                continue
            end
            prefactor, ms3 = ms_mult(MajoranaString(msum1.nfermions, ms1), MajoranaString(msum2.nfermions, ms2))
            add!(res, ms3, prefactor * tonumber(coeff1) * tonumber(coeff2))
        end
    end
    return res
end

function pop_id!(msum::MajoranaSum)
    if haskey(msum.Majoranas, 0)
        delete!(msum.Majoranas, 0)
    end
    return
end

function fock_filter(msum::MajoranaSum)
    clean_res = similar(msum)
    singles_filter = create_max_single_filter(nfermions(msum))
    for (ms, coeff) in msum.Majoranas
        if compute_max_single(ms, 2 * nfermions(msum), singles_filter) > 0
            continue
        end
        set!(clean_res, ms, coeff)
    end
    return clean_res
end

function overlap_with_fock(msum::MajoranaSum, fock_state; add_pref=0.)
    res = 0.
    singles_filter = create_max_single_filter(nfermions(msum))
    for (ms, coeff) in msum.Majoranas
        res += fockevaluate(ms, coeff, singles_filter, fock_state)
    end
    return res + add_pref
end

function fockevaluate(ms::TT, coeff, singles_filter, fock_state) where {TT<:Integer}
    if compute_max_single(ms, 2, singles_filter) > 0
        return 0.
    end
    num_pref = 0
    for site in fock_state
        num_pref += (ms >> (2 * site - 1)) & 1
    end
    ms_w = get_weight(ms)
    sign = (1im)^omega_L_mult(ms) * (1im)^(ms_w / 2) * (-1)^num_pref
    return tonumber(coeff) * sign
end

function overlap_with_fock_spinful(mslist, up_sites_with_particle, down_sites_with_particle, nsites; add_pref=0.)
    fock_state = []
    for up_site in up_sites_with_particle
        push!(fock_state, 2 * up_site - 1)
    end
    for down_site in down_sites_with_particle
        push!(fock_state, 2 * down_site)
    end
    return overlap_with_fock(mslist, fock_state; add_pref=add_pref)
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