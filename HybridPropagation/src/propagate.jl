
#Still not resolved why this is required
function PauliPropagation.PropagationBase.propagate!(circuit, prop_cache::AbstractHybridPropagationCache, params=nothing; kwargs...)
    circuit = reverse(circuit)
    if params isa Vector
        params = reverse(params)
    end

    return PauliPropagation.PropagationBase._propagate!(circuit, prop_cache, params; kwargs...)
end


function PauliPropagation.PropagationBase.propagate!(circ, hsum::HybridSum, thetas, fermions_filter, qubits_filter; max_weight=Inf, min_abs_coeff=1e-10, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...)
    # add filters to other kwargs as PauliPropagation.propagate! does not have those fields
    kwargs_dict::Dict{Symbol, Any} = Dict(kwargs)
    kwargs_dict[:qubits_filter] = qubits_filter
    kwargs_dict[:fermions_filter] = fermions_filter
    
    return PauliPropagation.propagate!(circ, hsum, thetas; max_weight=max_weight, min_abs_coeff=min_abs_coeff, max_freq=max_freq, max_sins=max_sins, customtruncfunc=customtruncfunc, kwargs_dict...)
end