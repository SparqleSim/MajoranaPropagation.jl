
function PropagationBase.applymergetruncate!(gate::FermionicRotation, prop_cache::MajoranaPropagationCache{MajoranaSum{TT, LookupCT}}, theta; truncate_each_mr=nothing, kwargs...) where {TT<:Integer}
    # get the Majorana strings and coefficients corresponding to the fermionic gate
    ms_rotations, coeffs, _ = MajoranaPropagation.getmajoranarotations(gate, nsites(prop_cache))

    # iterate over individual Majorana rotations and apply them to the Majorana sum
    for (gate_ms, coeff) in zip(ms_rotations, coeffs)
        # multiply coefficient by 2 since `::MajoranaRotation` implements exp(-i * theta/2 * mstring)
        applytoall!(gate_ms, prop_cache, theta * coeff * 2.0; kwargs...)

        # merge the auxiliary Majorana sum into the original one and empty the auxiliary one
        merge!(prop_cache; kwargs...)
    end

    return prop_cache
end

function _applycos(coeff::LookupCT{TT}, prefactor) where {TT<:Integer}
    new_coeff = deepcopy(coeff)
    for term in new_coeff.terms 
        push!(term.cos_theta_pref, prefactor)
    end
    return new_coeff
end
function _applysin(coeff::LookupCT{TT}, prefactor, gate_int, nfermions::Int) where {TT<:Integer}
    new_coeff = deepcopy(coeff)
    for term in new_coeff.terms 
        push!(term.sin_theta_pref, prefactor)
        pref, term.cumulative_string = ms_mult(gate_int, term.cumulative_string, nfermions)
        term.expression_pref *= 1im * pref 
    end
    return new_coeff
end


"""
The splitting rule for exp(i theta gate_string / 2) ms exp(-i theta gate_string / 2) is
-) ms, if [gate_string, ms] = 0
-) cos(theta) ms + i sin(theta) gate_string * ms, if {gate_string, ms} = 0
"""
function PropagationBase.applytoall!(gate::MajoranaRotation, prop_cache::MajoranaPropagationCache{MajoranaSum{TT, LookupCT}}, theta; kwargs...) where{TT<:Integer}
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    gate_int = gate.ms_int

    # loop over all Majorana strings and their coefficients in the Majorana sum
    for (ms_int, coeff) in msum
        if commutes(gate_int, ms_int)
            # if the gate commutes with the Majorana string, do nothing
            continue
        end

        # else we know the gate will split the Majorana string into two
        coeff1 = _applycos(coeff, theta)
        _, new_ms = ms_mult(gate_int, ms_int, nfermions(msum))
        coeff2 = _applysin(coeff, theta, gate_int, nfermions(msum))

        # set the coefficient of the original Majorana string
        set!(msum, ms_int, coeff1)

        # set the coefficient of the new Majorana string in the aux_psum
        # we can set the coefficient because MajoranaRotations create non-overlapping new Majorana strings
        set!(aux_msum, new_ms, coeff2)
    end

    return
end