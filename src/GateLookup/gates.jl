
function _evaluate_coeffs(coeff::LookupCT{TT}, ms_int::TT, nfermions::Int) where {TT<:Integer}
    res = 0.
    new_ms = TT(0)
    for term in coeff.terms
        pref, new_ms = ms_mult(term.cumulative_string, ms_int, nfermions)
        res += real(pref * term.expression_pref)
    end
    return res, new_ms
end


function PropagationBase.applymergetruncate!(gate::FermionicRotationLookup{TT}, prop_cache::MajoranaPropagationCache, theta; kwargs...) where {TT<:Integer}
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    # loop over all Majorana strings and their coefficients in the Majorana sum
    for (ms_int, coeff) in msum
        if !haskey(gate.lookup_table, ms_int)
            # if the gate commutes with the Majorana string, do nothing
            continue
        end

        propagated_strings = gate.lookup_table[ms_int]

        for (ms_key, gate_coeff) in propagated_strings
            if ms_key == ms_int
                all_cos_coeff, _ = _evaluate_coeffs(gate_coeff, ms_int, nfermions(msum))
                set!(msum, ms_int, all_cos_coeff)
            else
                resulting_coeff, new_string = _evaluate_coeffs(gate_coeff, ms_key, nfermions(msum))
                @show resulting_coeff, new_string
                add!(aux_msum, new_string, resulting_coeff)
            end
        end
    end

    return
end