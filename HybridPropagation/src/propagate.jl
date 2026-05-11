
#Still not resolved why this is required
function PauliPropagation.PropagationBase.propagate!(circuit, prop_cache::AbstractHybridPropagationCache, params=nothing; kwargs...)
    circuit = reverse(circuit)
    if params isa Vector
        params = reverse(params)
    end

    return PauliPropagation.PropagationBase._propagate!(circuit, prop_cache, params; kwargs...)
end
