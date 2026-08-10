"""
    two_body_term(nsites::Integer, indices, coeff::CT) where {CT}

Build the one-body (two fermion-operator) term of a quantum chemistry Hamiltonian for the orbital pair `indices` ``= (p, q)``, scaled by `coeff`: the spin-summed hopping ``\\sum_{\\sigma, \\sigma'} (f_{p\\sigma}^\\dagger f_{q\\sigma'} + f_{q\\sigma'}^\\dagger f_{p\\sigma})`` for ``p \\neq q``, or ``\\sum_{\\sigma, \\sigma'} f_{p\\sigma}^\\dagger f_{p\\sigma'}`` for ``p = q``.
"""
function two_body_term(nsites::Integer, indices, coeff::CT) where {CT}
    msum = MajoranaSum(CT, nsites, true)
    if indices[1] == indices[2]
        add!(msum, :nup, indices[1], coeff)
        add!(msum, :ndn, indices[1], coeff)
        add!(msum, :hop_on_site, indices[1], coeff)
    else
        add!(msum, :hopup, indices, coeff)
        add!(msum, :hopdn, indices, coeff)
        add!(msum, :hopupdn, indices, coeff)
        add!(msum, :hopupdn, reverse(indices), coeff)
    end
    return msum
end


"""
    four_body_term(nsites, p, q, r, s)

Build the two-body (four fermion-operator) term of a quantum chemistry Hamiltonian for the orbital indices ``(p, q, r, s)``: products ``f_{p}^\\dagger f_{q}^\\dagger f_{r} f_{s}`` with the spin patterns ``\\uparrow\\uparrow\\uparrow\\uparrow``, ``\\downarrow\\downarrow\\downarrow\\downarrow``, and ``\\uparrow\\downarrow\\uparrow\\downarrow``, summed over the 8 index permutations under which a real two-electron integral is invariant: ``h_{pqrs} = h_{prqs} = h_{sqrp} = h_{srqp} = h_{qpsr} = h_{qspr} = h_{rpsq} = h_{rspq}``.
Returns a `MajoranaSum` with real `Float64` coefficients.
"""
function four_body_term(nsites, p, q, r, s)
    op = MajoranaSum(ComplexF64, nsites, true)

    op += compute_permutatation_four_body_term(nsites, p, q, r, s)
    op += compute_permutatation_four_body_term(nsites, p, r, q, s)
    op += compute_permutatation_four_body_term(nsites, s, q, r, p)
    op += compute_permutatation_four_body_term(nsites, s, r, q, p)
    op += compute_permutatation_four_body_term(nsites, q, p, s, r)
    op += compute_permutatation_four_body_term(nsites, q, s, p, r)
    op += compute_permutatation_four_body_term(nsites, r, p, s, q)
    op += compute_permutatation_four_body_term(nsites, r, s, p, q)

    out_res = MajoranaSum(Float64, nsites, true)
    for (ms, coeff) in op
        if abs(coeff) != 0.0
            @assert abs(imag(coeff)) == 0.0
            set!(out_res, ms, real(coeff))
        end
    end
    return out_res
end

function compute_permutatation_four_body_term(nsites, p, q, r, s)
    op_up = MajoranaSum(ComplexF64, nsites, true)
    TT = getinttype(nfermions(op_up))
    add!(op_up, TT(0), 1.0)
    op_up *= MajoranaSum(nsites, :fupdag, p)
    op_up *= MajoranaSum(nsites, :fupdag, q)
    op_up *= MajoranaSum(nsites, :fup, r)
    op_up *= MajoranaSum(nsites, :fup, s)

    op_dn = MajoranaSum(ComplexF64, nsites, true)
    add!(op_dn, TT(0), 1.0)
    op_dn *= MajoranaSum(nsites, :fdndag, p)
    op_dn *= MajoranaSum(nsites, :fdndag, q)
    op_dn *= MajoranaSum(nsites, :fdn, r)
    op_dn *= MajoranaSum(nsites, :fdn, s)

    op_mixed = MajoranaSum(ComplexF64, nsites, true)
    add!(op_mixed, TT(0), 1.0)
    op_mixed *= MajoranaSum(nsites, :fupdag, p)
    op_mixed *= MajoranaSum(nsites, :fdndag, q)
    op_mixed *= MajoranaSum(nsites, :fup, r)
    op_mixed *= MajoranaSum(nsites, :fdn, s)
    full_op = MajoranaSum(ComplexF64, nsites, true)
    full_op += op_up
    full_op += op_dn
    full_op += op_mixed
    return full_op
end
