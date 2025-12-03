struct hybridRotation{TT<:Integer} <: PauliPropagation.ParametrizedGate
    hs::hybridString{TT}
end

struct HybridGate <: PauliPropagation.ParametrizedGate
    q_symbols::Vector{Symbol}
    q_sites::Vector{Int}
    f_symbol::Symbol
    f_sites::Vector{Int}
end

function gethybridrotations(gate::HybridGate, nqubits::Integer, nfermionic_sites::Integer)
    # construct hsum encoding the hybrid gate
    hsum = hybridSum(nqubits, gate.q_symbols, gate.q_sites, nfermionic_sites, gate.f_symbol, gate.f_sites)

    #remove coefficient associated to identity
    pop_id!(hsum)

    rotations::Vector{hybridRotation} = []
    coefficients::Vector{Float64} = []
    for (hs, coeff) in hsum
        push!(rotations, hybridRotation(hybridString(hsum.nqubits, nfermions(hsum), hs)))
        push!(coefficients, coeff)
    end

    return rotations, coefficients
end

function applytoall!(gate::hybridRotation, theta, hsum::hybridSum, aux_hsum::hybridSum; kwargs...)
    cos_val = cos(theta)
    sin_val = sin(theta)

    # separate gate into qubits and fermions part 
    Pj = gate.hs.string & kwargs[:qubits_filter]
    mu_v = gate.hs.string & kwargs[:fermions_filter]

    for (hs, coeff) in hsum  
        Pi = hs & kwargs[:qubits_filter]
        mu_w = hs & kwargs[:fermions_filter]

        if PauliPropagation._bitcommutes(Pi, Pj) == MajoranaPropagation.commutes(mu_v, mu_w)
            continue 
        else 
            coeff1 = _applycos(coeff, cos_val)

            sign, new_ms = ms_mult(mu_v, mu_w, nfermions(hsum))
            Pk, pk_sign = pauliprod(Pi, Pj)
            coeff2 = _applysin(coeff, sin_val * real(-1im * sign * pk_sign))
            hs2 = Pk | new_ms

            MajoranaPropagation.set!(hsum, hs, coeff1)
            MajoranaPropagation.set!(aux_hsum, hs2, coeff2)
        end
    end
    return 
end

function applymergetruncate!(gate::HybridGate, hsum::hybridSum{TT,CT}, aux_hsum::hybridSum{TT,CT}, thetas, param_idx; kwargs...) where {TT<:Integer,CT}
    # get the hybrid strings and coefficients corresponding to the hybrid gate
    hs_rotations, coeffs = gethybridrotations(gate, hsum.nqubits, hsum.nfermionic_sites)

    # get the current parameter
    theta = thetas[param_idx]

    # iterate over individual hybrid rotations and apply them to the hybrid sum
    for (gate_hs, coeff) in zip(hs_rotations, coeffs)
        # multiply coefficient by 2 since exponential implements exp(-i * theta/2 * hstring)
        applytoall!(gate_hs, theta * coeff * 2.0, hsum, aux_hsum; kwargs...)

        # merge the auxiliary Majorana sum into the original one and empty the auxiliary one
        hsum, aux_hsum = mergeandempty!(hsum, aux_hsum)

        # truncate after each Majorana rotation 
        checktruncationonall!(hsum; kwargs...)
    end

    # decrement parameter index during back propagation
    return hsum, aux_hsum, param_idx - 1
end

function PauliPropagation.propagate!(circ, hsum::hybridSum, thetas, qubits_filter, fermions_filter; max_weight=Inf, min_abs_coeff=1e-10, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...)
    # add filters to other kwargs as PauliPropagation.propagate! does not have those fields
    kwargs_dict::Dict{Symbol, Any} = Dict(kwargs)
    kwargs_dict[:qubits_filter] = qubits_filter
    kwargs_dict[:fermions_filter] = fermions_filter
    
    return PauliPropagation.propagate!(circ, hsum, thetas; max_weight=max_weight, min_abs_coeff=min_abs_coeff, max_freq=max_freq, max_sins=max_sins, customtruncfunc=customtruncfunc, kwargs_dict...)
end