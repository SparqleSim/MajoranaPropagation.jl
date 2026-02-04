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

function add!(ms::MajoranaSum{TT,CT}, symbol::Symbol, sites) where {TT<:Integer,CT}
    add!(ms, MajoranaSum(ms.nsites, symbol, sites))
    return ms
end

function add!(ms::MajoranaSum{TT,CT}, ms2::MajoranaSum{TT,CT}) where {TT<:Integer,CT}
    mergewith!(+, ms.Majoranas, ms2.Majoranas)
    return ms
end


function add!(ms::MajoranaSum{TT,CT}, ms2::MajoranaString{TT}, value::CT) where {TT<:Integer,CT}
    add!(ms, ms2.gammas, value)
end

function add!(ms::MajoranaSum, ms2_gammas::TT, value::CT) where {TT<:Integer,CT}
    if haskey(ms.Majoranas, ms2_gammas)
        ms.Majoranas[ms2_gammas] += value
    else
        ms.Majoranas[ms2_gammas] = value
    end
end

function set!(ms::MajoranaSum{TT,CT}, ms2::MajoranaString, value::CT) where {TT<:Integer,CT}
    set!(ms, ms2.gammas, value)
    return
end

function set!(ms::MajoranaSum{TT,CT}, ms2::TT, value::CT) where {TT<:Integer,CT}
    ms.Majoranas[ms2] = value
    return
end

function majoranas(ms::MajoranaSum)
    return keys(ms.Majoranas)
end

function coefficients(ms::MajoranaSum)
    return values(ms.Majoranas)
end

function nfermions(ms::MajoranaSum)
    if ms.is_spinful
        return 2 * ms.nsites
    else
        return ms.nsites
    end
end

function nfermions(ms::MajoranaString)
    return ms.nfermions
end

function Base.delete!(ms::MajoranaSum{TT,CT}, ms2::MajoranaString{TT}) where {TT<:Integer,CT}
    delete!(ms.Majoranas, ms2.gammas)
end
function Base.delete!(ms::MajoranaSum{TT,CT}, ms2_gammas::TT) where {TT<:Integer,CT}
    delete!(ms.Majoranas, ms2_gammas)
end

function Base.pop!(ms::MajoranaSum{TT,CT}, ms2_gammas::TT) where {TT<:Integer,CT}
    return pop!(ms.Majoranas, ms2_gammas, 0.)
end

function Base.length(ms::MajoranaSum)
    return length(ms.Majoranas)
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

function Base.show(io::IO, ms::MajoranaSum)
    max_display = 8
    print(io, "MajoranaSum with $(length(ms)) term(s):")
    for (i, (mstring, coeff)) in enumerate(ms.Majoranas)
        if i <= max_display
            print(io, "\n")
            print(io, "    $(coeff) * $(reverse(string(mstring; base=2, pad=2 * nfermions(ms))))")
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

function Base.:(*)(msum1::MajoranaSum{TT,CT1}, msum2::MajoranaSum{TT,CT2}) where {TT<:Integer,CT1,CT2}
    _checknfermions(msum1, msum2)
    res = MajoranaSum(ComplexF64, msum1.nsites, msum1.is_spinful)
    for (ms1, coeff1) in msum1
        for (ms2, coeff2) in msum2
            prefactor, ms3 = ms_mult(ms1, ms2, nfermions(msum1))
            add!(res, ms3, prefactor * coeff1 * coeff2)
        end
    end
    all_real = sum(abs.(imag.(coefficients(res)))) ≈ 0.
    #if all coefficients are real, convert back to real type and return that
    if all_real
        res_real = MajoranaSum(Float64, res.nsites, res.is_spinful)
        for (ms, coeff) in res
            set!(res_real, ms, real(coeff))
        end
        return res_real
    end

    return res
end

function Base.:(*)(coeff::CT, msum::MajoranaSum{TT,CT}) where {TT<:Integer,CT}
    res = similar(msum)
    for (ms1, coeff1) in msum
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

function commutator(msum1::MajoranaSum{TT,CT1}, msum2::MajoranaSum{TT,CT2}) where {TT<:Integer,CT1,CT2}
    res = MajoranaSum(ComplexF64, msum1.nsites, msum1.is_spinful)
    for (ms1, coeff1) in msum1
        for (ms2, coeff2) in msum2
            if commutes(ms1, ms2)
                continue
            end
            prefactor, ms3 = ms_mult(ms1, ms2, nfermions(msum1))
            add!(res, ms3, prefactor * coeff1 * coeff2)
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
    singles_filter = create_unpaired_mask(nfermions(msum))
    for (ms, coeff) in msum
        if compute_unpaired(ms, singles_filter) > 0
            continue
        end
        set!(clean_res, ms, coeff)
    end
    return clean_res
end

"""
    fockstate 
A struct to represent a Fock basis state. Contains a list of occupied sites and a boolean indicating whether the state is spinful or spinless.
"""
struct fockstate
    occupied_sites::Vector{Int}
    is_spinful::Bool
end

""" 
    fockstate(occupied_sites::Vector{Int})
Create a spinless Fock basis state given a list of occupied sites.
"""
function fockstate(occupied_sites::Vector)
    return fockstate(occupied_sites, false)
end

function fockstate(occupied_sites::StepRange)
    return fockstate(collect(occupied_sites))
end



""" 
    fockstate(up_occupied_sites::Vector{Int}, down_occupied_sites::Vector{Int})
Create a spinful Fock basis state given a list of occupied sites for spin-up and spin-down fermions.
"""
function fockstate(up_occupied_sites::Vector, down_occupied_sites::Vector)
    occupied_sites::Vector{Int} = []
    for site in up_occupied_sites
        push!(occupied_sites, 2 * site - 1)
    end
    for site in down_occupied_sites
        push!(occupied_sites, 2 * site)
    end
    sort!(occupied_sites)
    return fockstate(occupied_sites, true)
end

function fockstate(up_occupied_sites::StepRange, down_occupied_sites::StepRange)
    return fockstate(collect(up_occupied_sites), collect(down_occupied_sites))
end

""" 
    overlapwithfock(msum::MajoranaSum, fock_state::fockstate)
Compute the overlap <fock_state|msum|fock_state> where fock_state is a `fockstate` object.
"""
function overlapwithfock(msum::MajoranaSum, fock_state::fockstate)
    @assert msum.is_spinful == fock_state.is_spinful "The MajoranaSum and the fock_state must both be spinful or both spinless."
    res = 0.
    unpaired_mask = create_unpaired_mask(nfermions(msum))
    for (ms, coeff) in msum
        res += tonumber(coeff) * overlapwithfock(ms, unpaired_mask, fock_state)
    end
    return res
end

"""
    overlapwithfock(ms::TT, unpaired_mask::TT, fock_state::fockstate) where {TT<:Integer}
Compute the overlap <fock_state|ms|fock_state> where fock_state is a `fockstate` object.
"""
function overlapwithfock(ms::TT, unpaired_mask::TT, fock_state::fockstate) where {TT<:Integer}
    if compute_unpaired(ms, unpaired_mask) > 0
        return 0.
    end
    num_pref = 0
    for site in fock_state.occupied_sites
        num_pref += (ms >> (2 * site - 1)) & 1
    end
    ms_w = get_weight(ms)
    sign = (1im)^omega_L_mult(ms) * (1im)^(ms_w / 2) * (-1)^num_pref
    return real(sign)
end

""" 
    overlapwithfock(msum::MajoranaSum{TT,CT}, fock_state_1::fockstate, fock_state_2::fockstate) where {TT<:Integer,CT}

Evaluate the matrix element <fock_state_1|ms|fock_state_2> where fock_state_j are Fock basis states given as list of integers indicating which sites are occupied.
"""
function overlapwithfock(msum::MajoranaSum{TT,CT}, fock_state_1::fockstate, fock_state_2::fockstate) where {TT<:Integer,CT}
    @assert is_spinful(msum) == fock_state_1.is_spinful == fock_state_2.is_spinful "The MajoranaSum and the fock_states must both be spinful or both spinless."
    res = 0.
    n_fermions = nfermions(msum)
    for (ms, coeff) in msum
        res += coeff * overlapwithfock(ms, fock_state_1, fock_state_2, n_fermions)
        @show bitstring(ms), res
    end
    return res
end

""" 
    overlapwithfock(ms::TT, fock_state_1::fockstate, fock_state_2::fockstate) where {TT<:Integer}

Evaluate the matrix element <fock_state_1|ms|fock_state_2> where fock_state_j are Fock basis states given as list of integers indicating which sites are occupied.
"""
function overlapwithfock(ms::TT, fock_state_1::fockstate, fock_state_2::fockstate, n_fermions) where {TT<:Integer}
    res = (1im)^omega_L_mult(ms)
    for i = 1:n_fermions
        gamma = ((ms >> (2 * i - 2)) & TT(1))
        gamma_prime = ((ms >> (2 * i - 1)) & TT(1))
        if gamma == gamma_prime
            if (i in fock_state_2.occupied_sites) != (i in fock_state_1.occupied_sites)
                res *= 0.
                break
            else
                res *= (1im * (-1)^(i in fock_state_2.occupied_sites))^gamma
            end
        else
            if (i in fock_state_2.occupied_sites) == (i in fock_state_1.occupied_sites)
                res *= 0.
                break
            else
                res *= (1im * (-1)^(i in fock_state_2.occupied_sites))^gamma_prime * (-1)^(sum((j in fock_state_1.occupied_sites) for j = min(i + 1, n_fermions):n_fermions))
            end
        end
    end
    return res
end

""" 
    overlapwithfock(msum::MajoranaSum, sites_with_particle_superposition::Vector{fockstate}, superposition_coefficients::Vector{<:Union{Real,Complex}})
Compute the overlap <superposition|msum|superposition> where 
- superposition is given as a vector of Fock basis states `sites_with_particle_superposition`
- superposition_coefficients are the coefficients of the superposition (assumed normalized)
"""
function overlapwithfock(msum::MajoranaSum, sites_with_particle_superposition::Vector{fockstate}, superposition_coefficients::Vector{<:Union{Real,Complex}})
    # check normalization 
    @assert sum(abs2, superposition_coefficients) ≈ 1. "Superposition coefficients must be normalized."
    res = 0.
    unpaired_mask = create_unpaired_mask(nfermions(msum))

    for (ms, coeff) in msum
        for (sites_with_particle, superposition_coefficient) in zip(sites_with_particle_superposition, superposition_coefficients)
            res += coeff * abs(superposition_coefficient)^2 * overlapwithfock(ms, unpaired_mask, sites_with_particle)
        end

        for k1 = 1:length(sites_with_particle_superposition)
            for k2 = k1+1:length(sites_with_particle_superposition)
                fock1 = sites_with_particle_superposition[k1]
                fock2 = sites_with_particle_superposition[k2]
                superposition_coeff1 = superposition_coefficients[k1]
                superposition_coeff2 = superposition_coefficients[k2]
                res += 2. * real(coeff * conj(superposition_coeff1) * superposition_coeff2 * overlapwithfock(ms, fock1, fock2, nfermions(msum)))
            end
        end
    end
    return res
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