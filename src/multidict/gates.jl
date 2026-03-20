function PropagationBase.applymergetruncate!(gate::FermionicGate, prop_cache::MajoranaMultiPropagationCache, theta; truncate_each_mr=nothing, kwargs...)
    # get the Majorana strings and coefficients corresponding to the fermionic gate
    ms_rotations, coeffs, truncate_after_each_majrot = getmajoranarotations(gate, nsites(prop_cache))
    if !isnothing(truncate_each_mr)
        truncate_after_each_majrot = truncate_each_mr
    end

    # iterate over individual Majorana rotations and apply them to the Majorana sum
    for (gate_ms, coeff) in zip(ms_rotations, coeffs)

        # multiply coefficient by 2 since `::MajoranaRotation` implements exp(-i * theta/2 * mstring)
        applytoall!(gate_ms, prop_cache, theta * coeff * 2.0; kwargs...)

        # merge the auxiliary Majorana sum into the original one and empty the auxiliary one 
        merge!(prop_cache; kwargs...)

        # truncate after each Majorana rotation 
        if truncate_after_each_majrot
            truncate!(prop_cache; kwargs...)
        end
    end

    if !truncate_after_each_majrot
        truncate!(prop_cache; kwargs...)
    end

    return prop_cache
end

function PropagationBase.applytoall!(gate::MajoranaRotation, prop_cache::MajoranaMultiPropagationCache, theta; n_levels::Int, kwargs...)
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    gate_int = gate.ms_int

    msum_keys = collect(keys(msum))
    @threads for iw in eachindex(msum_keys)
    #for iw in eachindex(msum_keys)
        weight_key = msum_keys[iw]
        aux_entry = get(aux_msum, weight_key, nothing)
        if isnothing(aux_entry)
            weight_sector = _get_weight_from_key(weight_key)
            aux_entry = similar(msum, weight_sector, n_levels)
            aux_msum[weight_key] = aux_entry
        end
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
        try
            set!(aux_msum, dict_key, new_ms, coeff2)
        catch e
            @show weight_key
            @show dict_key
            @show bitstring(ms_int)
            @show bitstring(new_ms)
            println("Error occurred while setting aux_msum entry for key $dict_key")
            @show e 
            @show collect(keys(aux_msum.MultiMajoranas))
            @show compute_unpaired(new_ms, unpaired_mask)
            println(hgjffghjk)
        end
    end
    return
end

function _is_gaussian(gate_ms::MajoranaRotation{TT}) where {TT<:Integer}
    return get_weight(gate_ms.ms_int) == 2
end