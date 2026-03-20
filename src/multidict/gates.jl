function PropagationBase.applytoall!(gate::MajoranaRotation, prop_cache::MajoranaMultiPropagationCache, theta; n_levels::Int, kwargs...)
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    gate_int = gate.ms_int

    msum_keys = collect(keys(msum))
    # Create any missing aux entries up-front to avoid concurrent writes to `aux_msum`.
    for weight_key in msum_keys
        if !haskey(aux_msum, weight_key)
            weight_sector = _get_weight_from_key(weight_key)
            aux_msum[weight_key] = similar(msum, weight_sector, n_levels)
        end
    end

    @threads for iw in eachindex(msum_keys)
        weight_key = msum_keys[iw]
        aux_entry = aux_msum[weight_key]
        _applymajoranarotation!(
            msum.MultiMajoranas[weight_key],
            aux_entry,
            gate_int,
            theta,
            nfermions(msum);
            weight_key,
            kwargs...,
        )
    end
    return prop_cache
end

function _applymajoranarotation!(msum_dict::Dict{TT,CT}, aux_msum::MajoranaSumMulti, gate_int, theta, n_fermions; level_mapper::Function,weight_key,unpaired_mask, kwargs...) where {TT<:Integer,CT}
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

        # set the coefficient of the new Pauli string in the corresponding aux_psum
        # we can set the coefficient because PauliRotations create non-overlapping new Pauli strings
        weight = get_weight(new_ms)
        dict_key = "$weight-$(level_mapper(new_ms))"
        set!(aux_msum, dict_key, new_ms, coeff2)
    end
    return
end

function _is_gaussian(gate_ms::MajoranaRotation{TT}) where {TT<:Integer}
    return get_weight(gate_ms.ms_int) == 2
end