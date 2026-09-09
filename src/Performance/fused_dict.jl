###
##
# Fused `applymergetruncate!` for the dictionary-based `MajoranaSum` cache, walking and mutating
# the backing `Dict` through its internal slot array as `PauliPropagation.Performance` does for
# `PauliSum`. The internals used here are the ones `PauliPropagation.Performance._check_dict_internals`
# verifies when PauliPropagation loads.
# A term whose coefficient falls below `min_abs_coeff` after the cosine scaling may still receive
# the sine contribution of its partner string, so it is only removed in a third pass over just
# those keys once every product has landed. The results are identical to the default path.
##
###

"""
    applymergetruncate!(gate::MajoranaRotation, prop_cache::MajoranaPropagationCache{<:MajoranaSum}, theta; fused::Bool=false, kwargs...)

Fused overload that truncates during gate application by walking and mutating the backing `Dict` through its internal slot array -- see file header.
`force_truncate=true` tests every term against the truncations, including those the gate leaves alone.
"""
function PropagationBase.applymergetruncate!(gate::MajoranaRotation, prop_cache::MajoranaPropagationCache{<:MajoranaSum}, theta;
    fused::Bool=false,
    min_abs_coeff::Real=1e-10, max_weight::Real=Inf, max_unpaired::Real=Inf, unpaired_mask=nothing,
    max_freq::Real=Inf, max_sins::Real=Inf, customtruncfunc=nothing,
    force_truncate::Bool=false, kwargs...)

    if !fused
        return invoke(PropagationBase.applymergetruncate!,
            Tuple{MajoranaRotation,MP.AbstractMajoranaPropagationCache,typeof(theta)},
            gate, prop_cache, theta;
            min_abs_coeff, max_weight, max_unpaired, unpaired_mask, max_freq, max_sins, customtruncfunc, kwargs...)
    end

    msum = storage(mainsum(prop_cache))
    gate_ms = gate.ms_int
    gate_ms_ps = MP.compute_parity_bits_and_shift(gate_ms, 2 * nfermions(prop_cache))
    omega_l_gate = MP.omega_L_mult(gate_ms)
    cos_val, sin_val = cos(theta), sin(theta)

    mask = MP._unpairedmask(prop_cache, unpaired_mask, max_unpaired)
    truncfunc = MP._truncationfunc(mask; min_abs_coeff, max_weight, max_unpaired, max_freq, max_sins, customtruncfunc)

    touched = Tuple{keytype(msum),valtype(msum)}[]
    # bounded by one push per filled, branching slot
    sizehint!(touched, length(msum))
    # the scaled terms that fail a truncation on their own; their partner's product may still
    # merge into them, so they are judged only after pass 2
    tentative = keytype(msum)[]

    # pass 1: walk the slots, scaling the branching terms in place by slot index
    nslots = length(msum.slots)
    @inbounds for i in 1:nslots
        Base.isslotfilled(msum, i) || continue

        term = msum.keys[i]
        if MP._commutes_evengate(term, gate_ms)
            # a product never lands on a commuting term, so this test is final
            force_truncate && truncfunc(term, msum.vals[i]) && Base._delete!(msum, i)
            continue
        end
        coeff = msum.vals[i]

        coeff1 = MP._applycos(coeff, cos_val)
        msum.vals[i] = coeff1
        truncfunc(term, coeff1) && push!(tentative, term)

        new_term, sign = MP._rotationproduct_evengate(term, gate_ms, gate_ms_ps, omega_l_gate)
        MP._truncateterm(new_term, max_weight, max_unpaired, mask) && continue
        push!(touched, (new_term, MP._applysin(coeff, sin_val * sign)))
    end

    # pass 2: fold each product into its final value and truncate it right away. `touched` holds
    # no duplicate keys, since XOR with a fixed gate string is a bijection on Majorana strings, so
    # every touched key gets exactly one, complete update here.
    for (new_term, delta) in touched
        index, sh = Base.ht_keyindex2_shorthash!(msum, new_term)
        if index > 0
            new_val = mergefunc(msum.vals[index], delta)
            if truncfunc(new_term, new_val)
                Base._delete!(msum, index)
            else
                msum.vals[index] = new_val
            end
        elseif !truncfunc(new_term, delta)
            Base._setindex!(msum, delta, new_term, -index, sh)
        end
    end

    # pass 3: the tentative terms that pass 2 did not settle (a term it updated was judged on its
    # merged value there, and comes out the same way here)
    for term in tentative
        index, _ = Base.ht_keyindex2_shorthash!(msum, term)
        index > 0 && truncfunc(term, msum.vals[index]) && Base._delete!(msum, index)
    end

    return prop_cache
end
