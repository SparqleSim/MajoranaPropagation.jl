#=function applymergetruncate!(gate::FermionicGate, msum::MajoranaSumMulti{TT, CT}, dd, thetas, param_idx; max_weight=Inf, min_abs_coeff=1e-10, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...) where {TT<:Integer,CT}
    # Pick out the next theta if gate is a ParametrizedGate.
    # Else set the paramter to nothing for clarity that theta is not used.
    if gate isa ParametrizedGate
        theta = thetas[param_idx]
        # If the gate is parametrized, decrement theta index by one.
        param_idx -= 1
    else
        theta = nothing
    end

    ms_rotations, coeffs = getmajoranarotations(gate, msum.nsites)

    # Apply the gate to all Pauli strings in psum, potentially writing into auxillary aux_psum in the process.
    # The pauli sums will be changed in-place

    #=dicts_lengths = [length(msum.MultiMajoranas[k]) for k in msum_keys]
    sorting_indices = sortperm(dicts_lengths; rev=true)

    to_submit = []
    merge_in_apply=true
    for idx in sorting_indices
        push!(to_submit, [gate, theta, msum.MultiMajoranas[msum_keys[idx]], all_aux_msums[msum_keys[idx]], merge_in_apply, msum_keys[idx]])
    end

    function submit(iter_all)
        applytoall!(iter_all[1:end-2]...; merge_sector=iter_all[end-1], weight_key=iter_all[end], kwargs...)
    end

    #ThreadPools.qbforeach(iter -> submit(iter), to_submit)
    ThreadPools.qforeach(iter -> submit(iter), to_submit)=#

    merge_in_apply=true
    for (gate_ms, coeff) in zip(ms_rotations, coeffs)
        all_aux_msums = Dict()
        for weight_key in keys(msum.MultiMajoranas)
            all_aux_msums[weight_key] = similar(msum, weight_key)
        end
        msum_keys = collect(keys(msum.MultiMajoranas))

        @threads for iw=1:length(msum_keys)
            weight_key = msum_keys[iw]
            aux_psum = all_aux_msums[weight_key]
            # multiply coefficient by 2 since exponential implements exp(-i * theta/2 * mstring)
            #@show typeof(gate_ms), typeof(msum.MultiMajoranas[weight_key]), typeof(aux_psum)
            applytoall!(gate_ms, theta * coeff * 2., msum.MultiMajoranas[weight_key], aux_psum; merge_sector=merge_in_apply, weight_key=weight_key, kwargs...)
            #println("----")
            #println(gyfhj)
            #@show weight_key
            #@show msum_weight
            #@show aux_psum
        end
        
        msum, all_aux_msums = mergeandempty!(msum, all_aux_msums; merge_sector=!merge_in_apply)
        for weight_key in keys(msum.MultiMajoranas)
            checktruncationonall!(MajoranaSum(msum.nsites, msum.is_spinful, msum.MultiMajoranas[weight_key]); max_weight, min_abs_coeff, max_freq, max_sins, customtruncfunc)
        end
    end
    return msum, dd, param_idx
end=#

function PropagationBase.applytoall!(gate::MajoranaRotation, prop_cache::MajoranaMultiPropagationCache, theta; merge_sector=false, kwargs...)
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    gate_int = gate.ms_int

    msum_keys = collect(keys(msum))
    @threads for iw in eachindex(msum_keys)
    #for iw in eachindex(msum_keys)
        weight_key = msum_keys[iw]
        aux_entry = get(aux_msum, weight_key, nothing)
        if isnothing(aux_entry)
            aux_entry = similar(msum, weight_key)
            aux_msum[weight_key] = aux_entry
        end
        _applymajoranarotation!(
            msum.MultiMajoranas[weight_key],
            aux_entry,
            gate_int,
            theta,
            nfermions(msum);
            merge_sector=merge_sector,
            kwargs...,
        )
    end
    return prop_cache
end

function _applymajoranarotation!(msum_dict::Dict{TT,CT}, aux_msum::MajoranaSumMulti, gate_int, theta, n_fermions; merge_sector=false, kwargs...) where {TT<:Integer,CT}
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
        set!(aux_msum, weight, new_ms, coeff2)
    end

    if merge_sector
        # merge aux_msum back into msum
        mergewith!(+, msum, aux_msum.MultiMajoranas[kwargs[:weight_key]])
        # empty aux_msum
        empty!(aux_msum.MultiMajoranas[kwargs[:weight_key]])
    end

    return
end