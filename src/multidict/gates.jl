using TimerOutputs

function PropagationBase.applymergetruncate!(gate::FermionicGate, prop_cache::MajoranaMultiPropagationCache, theta; to, truncate_each_mr=nothing, kwargs...)
    # get the Majorana strings and coefficients corresponding to the fermionic gate
    ms_rotations, coeffs, truncate_after_each_majrot = getmajoranarotations(gate, nsites(prop_cache))
    if !isnothing(truncate_each_mr)
        truncate_after_each_majrot = truncate_each_mr
    end

    # iterate over individual Majorana rotations and apply them to the Majorana sum
    for (gate_ms, coeff) in zip(ms_rotations, coeffs)

        # multiply coefficient by 2 since `::MajoranaRotation` implements exp(-i * theta/2 * mstring)
        @timeit to "applytoall" applytoall!(gate_ms, prop_cache, theta * coeff * 2.0; kwargs...)

        # merge the auxiliary Majorana sum into the original one and empty the auxiliary one 
        @timeit to "merge" merge!(prop_cache; to, kwargs...)

        # truncate after each Majorana rotation 
        if truncate_after_each_majrot
            @timeit to "truncate" truncate!(prop_cache; kwargs...)
        end
    end

    if !truncate_after_each_majrot
        @timeit to "truncate" truncate!(prop_cache; kwargs...)
    end

    return prop_cache
end

function PropagationBase.applytoall!(gate::MajoranaRotation, prop_cache::MajoranaMultiPropagationCache, theta; kwargs...)
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    gate_int = gate.ms_int

    @threads for iw in eachindex(msum)
        _applymajoranarotation!(
            msum.MultiMajoranas[iw],
            aux_msum.MultiMajoranas[iw],
            gate_int,
            theta,
            nfermions(msum);
            kwargs...,
        )
    end
    return prop_cache
end

function _applymajoranarotation!(msum_dict::Dict{TT,CT}, aux_dict::Dict{TT,CT}, gate_int, theta, n_fermions; kwargs...) where {TT<:Integer,CT}
    cos_val = cos(theta)
    sin_val = sin(theta)
    
    # loop over all Pauli strings and their coefficients in the Pauli sum
    for (ms_int, coeff) in msum_dict
        if commutes(gate_int, ms_int)
            # if the gate commutes with the pauli string, do nothing
            continue
        end

        # else we know the gate will split th Pauli string into two
        coeff1 = _applycos(coeff, cos_val)
        sign, new_ms = ms_mult(gate_int, ms_int, n_fermions)
        coeff2 = _applysin(coeff, sin_val * real((-1im) * sign))

        # set the coefficient of the original Pauli string
        msum_dict[ms_int] = coeff1
        # add the new Pauli string with its coefficient
        aux_dict[new_ms] = coeff2
    end
    return
end

function _is_gaussian(gate_ms::MajoranaRotation{TT}) where {TT<:Integer}
    return get_weight(gate_ms.ms_int) == 2
end