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
    singles_mask = create_max_single_filter(nfermions(msum))
    for (ms, coeff) in msum.Majoranas
        res += fockevaluate(ms, coeff, singles_mask, fock_state)
    end
    return res + add_pref
end

function fockevaluate(ms::TT, coeff, singles_mask, fock_state) where {TT<:Integer}
    if compute_max_single(ms, 2, singles_mask) > 0
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