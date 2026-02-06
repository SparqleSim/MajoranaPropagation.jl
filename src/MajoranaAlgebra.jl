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

function fock_mask(msum::MajoranaSum)
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
A struct to represent a Fock basis state.
"""
struct fockstate{TT<:Integer}
    n_sites::Int
    is_spinful::Bool
    occupied_sites::TT
end

"""
    fockstate(n_sites::Int, occupied_sites_list::Vector{Int})
Create a spinless Fock basis state given a list of occupied sites.
"""
function fockstate(n_sites::Int, occupied_sites_list::Vector)
    TT = getinttype(n_sites)
    occupied_sites = _bitonesat(TT, (2 * site -1 for site in occupied_sites_list))
    return fockstate(n_sites, false, occupied_sites)
end

"""
    fockstate(n_sites::Int, up_occupied_sites_list::Vector, down_occupied_sites_list::Vector)
Create a spinful Fock basis state given a list of occupied sites for spin-up and spin-down fermions.
"""
function fockstate(n_sites::Int, up_occupied_sites_list::Vector, down_occupied_sites_list::Vector)
    TT = getinttype(2 * n_sites)
    occupied_sites_list::Vector{Int} = []
    for site in up_occupied_sites_list
        push!(occupied_sites_list, 2 * site - 1)
    end
    for site in down_occupied_sites_list
        push!(occupied_sites_list, 2 * site)
    end
    sort!(occupied_sites_list)
    occupied_sites = _bitonesat(TT, (2 * site -1 for site in occupied_sites_list))
    return fockstate(n_sites, true, occupied_sites)
end

"""
    overlapwithfock(msum::MajoranaSum, fock_state::fockstate)
Compute the overlap <fock_state|msum|fock_state> where fock_state is a `fockstate` object.
"""
function overlapwithfock(msum::AbstractMajoranaSum, fock_state::fockstate)
    @assert is_spinful(msum) == fock_state.is_spinful "The MajoranaSum and the fock_state must both be spinful or both spinless."
    res = 0.
    unpaired_mask = create_unpaired_mask(nfermions(msum))
    for (ms, coeff) in msum
        res += tonumber(coeff) * overlapwithfock(ms, unpaired_mask, fock_state)
    end
    return res
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
    overlapwithfock(ms::TT, unpaired_mask::TT, fock_state::fockstate) where {TT<:Integer}
Compute the overlap <fock_state|ms|fock_state> where fock_state is a `fockstate` object.
"""
function overlapwithfock(ms::TT, unpaired_mask::TT, fock_state::fockstate) where {TT<:Integer}
    if compute_unpaired(ms, unpaired_mask) > 0
        return 0.
    end
    number_pref = get_weight(ms & fock_state.occupied_sites)
    sign = (1im)^omega_L_mult(ms) * (1im)^(get_weight(ms) / 2) * (-1)^number_pref
    return real(sign)
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
# TODO: fix, now broken 
function overlapwithfock(msum::MajoranaSum, sites_with_particle_superposition::Vector{fockstate}, superposition_coefficients::Vector{<:Union{Real,Complex}})
    @error "Currently non supported"
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