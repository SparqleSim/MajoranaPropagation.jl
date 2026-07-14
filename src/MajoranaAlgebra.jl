# TODO: check if this makes sense
const _parity_mask_cache = Dict{Tuple{DataType, Int}, Any}()

function _get_parity_masks(::Type{TT}, Nbits::Int) where {TT<:Integer}
    key = (TT, Nbits)
    cached = get(_parity_mask_cache, key, nothing)
    if cached === nothing
        full_mask = (TT(1) << Nbits) - TT(1)
        half_mask = full_mask >> 1
        cached = (full_mask, half_mask)
        _parity_mask_cache[key] = cached
    end
    return cached::Tuple{TT, TT}
end

function compute_parity_bits_and_shift(u::TT, Nbits::Int) where {TT<:Integer}
    if Nbits <= 1
        return TT(0)
    end

    full_mask, half_mask = _get_parity_masks(TT, Nbits)

    p = u & half_mask

    s = 1
    while s < Nbits
        p ⊻= (p << s)
        s <<= 1
    end

    p = p << 1
    return p & full_mask
end

function omega_L_mult(ms1::MajoranaString, ms2::MajoranaString)
    return omega_L_mult_right_string_parity_shifted(ms1.gammas, compute_parity_bits_and_shift(ms2.gammas, 2 * ms1.nfermions))
end

function omega_L_mult(ms1::TT, ms2::TT, Nbits) where {TT<:Integer}
    return omega_L_mult_right_string_parity_shifted(ms1, compute_parity_bits_and_shift(ms2, Nbits))
end

function omega_L_mult_right_string_parity_shifted(ms1::TT, ms2_ps::TT) where {TT<:Integer}
    return mod(get_weight(ms1 & ms2_ps), 2)
end

function omega_L_mult(ms::TT) where {TT<:Integer}
    return (get_weight(ms) >> 1) & 1
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

function majoranarotationproduct(term_ms::TT, gate_ms::TT, gate_ms_ps::TT) where {TT<:Integer}
    new_term = gate_ms ⊻ term_ms
    omega_l_gate = omega_L_mult(gate_ms)
    omega_l_term = omega_L_mult(term_ms)
    omega_term_gate = omega_mult(term_ms, gate_ms)
    omega_l_term_gate = omega_L_mult_right_string_parity_shifted(term_ms, gate_ms_ps)
    fprefactor_precomputed = omega_l_gate * omega_l_term + omega_term_gate * (omega_l_gate + omega_l_term + 1)
    exp = omega_l_term_gate + fprefactor_precomputed
    sign = iseven(exp) ? 1 : -1
    return new_term, sign
end

"""
    _commutes_evengate(term_ms, gate_ms)

Commutation check specialized to even-weight `gate_ms`: the weight product in `omega_mult` is then
even for any term, so commutation reduces to an even popcount of the overlap.
"""
@inline function _commutes_evengate(term_ms::TT, gate_ms::TT) where {TT<:Integer}
    return iseven(count_ones(term_ms & gate_ms))
end

"""
    _rotationproduct_evengate(term_ms, gate_ms, gate_ms_ps, omega_l_gate)

`majoranarotationproduct` specialized to even-weight `gate_ms` and anticommuting `term_ms`
"""
@inline function _rotationproduct_evengate(term_ms::TT, gate_ms::TT, gate_ms_ps::TT, omega_l_gate::Int) where {TT<:Integer}
    new_term = gate_ms ⊻ term_ms
    omega_l_term = (count_ones(term_ms) >> 1) & 1
    omega_l_term_gate = count_ones(term_ms & gate_ms_ps) & 1
    exp = omega_l_term_gate + omega_l_gate * omega_l_term + omega_l_gate + omega_l_term + 1
    sign = iseven(exp) ? 1 : -1
    return new_term, sign
end

"""
    _rotationproduct_evengate_commuting(term_ms, gate_ms, gate_ms_ps, omega_l_gate)

Sign and result of the product `gate_ms * term_ms` specialized to even-weight `gate_ms` and commuting `term_ms` (the imaginary-time splitting branch)
"""
@inline function _rotationproduct_evengate_commuting(term_ms::TT, gate_ms::TT, gate_ms_ps::TT, omega_l_gate::Int) where {TT<:Integer}
    new_term = gate_ms ⊻ term_ms
    omega_l_term = (count_ones(term_ms) >> 1) & 1
    omega_l_term_gate = count_ones(term_ms & gate_ms_ps) & 1
    exp = omega_l_term_gate + omega_l_gate * omega_l_term
    sign = iseven(exp) ? 1 : -1
    return new_term, sign
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