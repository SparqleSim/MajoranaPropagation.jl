struct MajoranaSumMulti{TT<:Integer,CT}
    nsites::Int
    is_spinful::Bool
    MultiMajoranas::Dict{Int, Dict{TT,CT}}
end

function MajoranaSumMulti(msum::MajoranaSum{TT,CT}) where {TT<:Integer,CT}
    nsites = msum.nsites
    is_spinful = msum.is_spinful
    multimajs = Dict{Int, Dict{TT,CT}}()

    for (ms_int, coeff) in msum.Majoranas
        ms_weight = get_weight(ms_int)
        if !haskey(multimajs, ms_weight)
            multimajs[ms_weight] = Dict{TT, CT}()
        end
        multimajs[ms_weight][ms_int] = coeff
    end
    return MajoranaSumMulti{TT,CT}(nsites, is_spinful, multimajs)
end

function similar(msum::MajoranaSumMulti{TT,CT}, W) where {TT<:Integer,CT}
    return MajoranaSumMulti(msum.nsites, msum.is_spinful,
                            Dict(W-2  => Dict{TT,CT}(), 
                                W => Dict{TT,CT}(),
                                W+2 => Dict{TT,CT}()))
end

function coefftype(msum::MajoranaSumMulti{TT,CT}) where {TT,CT}
    return CT
end

function Base.length(msum::MajoranaSumMulti{TT,CT}) where {TT<:Integer,CT}
    total_strings = 0
    for (weight_key, dict) in msum.MultiMajoranas
        nstrings = length(dict)
        total_strings += nstrings
    end
    return total_strings
end

function show_stats(msum::MajoranaSumMulti{TT,CT}) where {TT<:Integer,CT}
    total_strings = 0
    for (weight_key, dict) in msum.MultiMajoranas
        nstrings = length(dict)
        total_strings += nstrings
    end
    for (weight_key, dict) in msum.MultiMajoranas
        nstrings = length(dict)
        println("Weight $weight_key: $nstrings strings ($(round(100. * nstrings / total_strings))%)")
    end
    println("Total strings: $total_strings")
end

function nfermions(ms::MajoranaSumMulti)
    if ms.is_spinful
        return 2 * ms.nsites
    else
        return ms.nsites
    end
end


function propagate!(circ, msum::MajoranaSumMulti{TT, CT}, thetas=nothing; max_weight=Inf, min_abs_coeff=1e-10, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...) where {TT<:Integer,CT}

    # check that max_freq and max_sins are only used a PathProperties type tracking them
    _checkfreqandsinfields(msum, max_freq, max_sins)

    # if circ is actually a single gate, promote it to a list [gate]
    # similarly the theta if it is a single number
    circ, thetas = _promotecircandthetas(circ, thetas)

    # if thetas is nothing, the circuit must contain only StaticGates
    # also check if the length of thetas equals the number of parametrized gates
    _checkcircandthetas(circ, thetas)

    # start from the last parameter if thetas is not nothing
    param_idx = thetas === nothing ? nothing : length(thetas)

    # get our auxillary Pauli sum that we will move splitting Pauli strings into 
    aux_psum = []

    ## TODO:
    # - decide where to reverse the circuit
    # - verbose option  
    # - more elegant param_idx incrementation
    for gate in reverse(circ)
        #@show gate, typeof(gate)
        #@show typeof(msum)
        msum, aux_psum, param_idx = applymergetruncate!(gate, msum, aux_psum, thetas, param_idx; max_weight, min_abs_coeff, max_freq, max_sins, customtruncfunc, kwargs...)
    end
    # TODO: potential bug: If there are Clifford gates in the circuit, merging may swap the psums.
    #                      Thi smeans that the original psum is not the one that is returned, and that the original psum is empty.
    return msum
end


function applymergetruncate!(gate::FermionicGate, msum::MajoranaSumMulti{TT, CT}, dd, thetas, param_idx; max_weight=Inf, min_abs_coeff=1e-10, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...) where {TT<:Integer,CT}
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
end


function applytoall!(gate::MajoranaRotation{TT}, theta, msum::Dict{TT,CT}, aux_msum::MajoranaSumMulti{TT,CT}; merge_sector=false, kwargs...) where {TT<:Integer,CT}
    cos_val = cos(theta)
    sin_val = sin(theta)

    gate_int = gate.ms.gammas

    # loop over all Pauli strings and their coefficients in the Pauli sum
    for (ms_int, coeff) in msum
        if commutes(gate_int, ms_int)
            # if the gate commutes with the pauli string, do nothing
            continue
        end

        # else we know the gate will split th Pauli string into two
        coeff1 = _applycos(coeff, cos_val)
        sign, new_ms = ms_mult(gate_int, ms_int, nfermions(aux_msum))
        coeff2 = _applysin(coeff, sin_val * real((-1im) * sign))

        # set the coefficient of the original Pauli string
        #set!(msum, ms_int, coeff1)
        msum[ms_int] = coeff1

        # set the coefficient of the new Pauli string in the corresponding aux_psum
        # we can set the coefficient because PauliRotations create non-overlapping new Pauli strings
        weight = get_weight(new_ms)
        set!(aux_msum.MultiMajoranas[weight], new_ms, coeff2)
    end

    if merge_sector
        # merge aux_msum back into msum
        mergewith!(+, msum, aux_msum.MultiMajoranas[kwargs[:weight_key]])
        # empty aux_msum
        empty!(aux_msum.MultiMajoranas[kwargs[:weight_key]])
    end

    return
end

function mergeandempty!(msums::MajoranaSumMulti{TT,CT}, aux_psum; merge_sector=true) where {TT<:Integer,CT}
    #@show "hi"
    sorted_keys = sort(collect(keys(msums.MultiMajoranas)))
    #@show sorted_keys

    #do the Delta W = 0 merges first
    if merge_sector 
        for weight_key in sorted_keys
            mergewith!(+,msums.MultiMajoranas[weight_key], aux_psum[weight_key].MultiMajoranas[weight_key])
            empty!(aux_psum[weight_key].MultiMajoranas[weight_key])
        end
    end

    # do the Delta W = +2 merges 
    #@threads for weight_key in sorted_keys
    for weight_key in sorted_keys
        if haskey(aux_psum, weight_key) == false
            @show sorted_keys 
            println(gfhj)
        end
        #@show weight_key
        if haskey(aux_psum[weight_key].MultiMajoranas, weight_key+2) == false
            @show sorted_keys 
            println(gfhjdfghjk)
        end
        if length(aux_psum[weight_key].MultiMajoranas[weight_key+2])>0
            if haskey(msums.MultiMajoranas, weight_key+2) == false
                #create new dict in msusm 
                msums.MultiMajoranas[weight_key+2] = Dict{TT,CT}()
            end
            mergewith!(+, msums.MultiMajoranas[weight_key+2], aux_psum[weight_key].MultiMajoranas[weight_key+2])
            empty!(aux_psum[weight_key].MultiMajoranas[weight_key+2])
        end
    end

    # do the Delta W = -2 merges
    @threads for weight_key in sorted_keys
        if length(aux_psum[weight_key].MultiMajoranas[weight_key-2])>0
            if haskey(msums.MultiMajoranas, weight_key-2) == false
                #create new dict in msusm 
                msums.MultiMajoranas[weight_key-2] = Dict{TT,CT}()
            end
            mergewith!(+, msums.MultiMajoranas[weight_key-2], aux_psum[weight_key].MultiMajoranas[weight_key-2])
            empty!(aux_psum[weight_key].MultiMajoranas[weight_key-2])
        end
    end

    sorted_keys_final = sort(collect(keys(msums.MultiMajoranas)))
    #remove empty dicts
    for weight_key in sorted_keys_final
        if length(msums.MultiMajoranas[weight_key])==0
            delete!(msums.MultiMajoranas, weight_key)
        end
    end

    return msums, aux_psum
end