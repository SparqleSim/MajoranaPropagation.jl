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

function Base.:(*)(msum1::MajoranaSum{TT,CT1}, msum2::MajoranaSum{TT,CT2}) where {TT<:Integer,CT1,CT2}
    _checknfermions(msum1, msum2)
    res = MajoranaSum(ComplexF64, msum1.nsites, msum1.is_spinful)
    for (ms1, coeff1) in zip(terms(msum1), coefficients(msum1))
        for (ms2, coeff2) in zip(terms(msum2), coefficients(msum2))
            prefactor, ms3 = ms_mult(ms1, ms2, nfermions(msum1))
            add!(res, ms3, prefactor * coeff1 * coeff2)
        end
    end
    all_real = maximum(abs.(imag.(coefficients(res)))) ≈ 0.
    #if all coefficients are real, convert back to real type and return that
    if all_real
        res_real = MajoranaSum(Float64, res.nsites, res.is_spinful)
        for (ms, coeff) in zip(terms(res), coefficients(res))
            set!(res_real, ms, real(coeff))
        end
        return res_real
    end
    return res
end

function Base.:(*)(coeff::Number, msum::MajoranaSum{TT,CT}) where {TT<:Integer,CT}
    res = similar(msum)
    for (ms1, coeff1) in zip(terms(msum), coefficients(msum))
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
    res = MajoranaSum(ComplexF64, nsites(msum1), is_spinful(msum1))
    for (ms1, coeff1) in zip(terms(msum1), coefficients(msum1))
        for (ms2, coeff2) in zip(terms(msum2), coefficients(msum2))
            if commutes(ms1, ms2)
                continue
            end
            prefactor, ms3 = ms_mult(ms1, ms2, nfermions(msum1))
            add!(res, ms3, prefactor * coeff1 * coeff2)
        end
    end
    return res
end

function scalarproduct(msum1::AbstractMajoranaSum, msum2::AbstractMajoranaSum)
    res = zero(eltype(coefficients(msum1)))
    if length(msum1) > length(msum2)
        # if msum1 has more terms than msum2, it's more efficient to loop over msum2 and check for each term if it is in msum1
        return scalarproduct(msum2, msum1)
    end
    for (ms1, coeff1) in zip(terms(msum1), coefficients(msum1))
        coeff2 = getmergedcoeff(msum2, ms1)
        res += coeff1 * coeff2
    end
    return res
end