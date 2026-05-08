#=
Not used in any capacity:
    max_weight (no definition of global weight)
    max_freq
    max_sins
=#

function PauliPropagation.PropagationBase.truncate!(
    prop_cache::AbstractHybridPropagationCache;
    max_weight::Real=Inf, min_abs_coeff=1e-10,
    max_pauli_weight::Real=Inf, max_majorana_weight::Real=Inf, max_unpaired::Real=Inf,
    max_freq::Real=Inf, max_sins::Real=Inf,
    unpaired_mask=nothing,
    customtruncfunc=nothing,
    fermion_filter=nothing, qubit_filter=nothing, 
    kwargs...

)
    if isnothing(unpaired_mask)
        unpaired_mask = create_unpaired_mask(nfermions(mainsum(prop_cache)))
    end
    if isnothing(fermion_filter) || isnothing(qubit_filter)
        fermion_filter, qubit_filter = create_filters(mainsum(prop_cache))
    end 

    function truncfunc(hstr, coeff)
        is_truncated = false
        if PauliPropagation.truncatemincoeff(coeff, min_abs_coeff)
            is_truncated = true
        elseif MajoranaPropagation.truncateunpaired(typeof(unpaired_mask)(hstr&fermion_filter), max_unpaired, unpaired_mask)
            is_truncated = true
        elseif MajoranaPropagation.truncatemajoranaweight((hstr&fermion_filter), max_majorana_weight)
            is_truncated = true
        elseif PauliPropagation.truncateweight((hstr&qubit_filter), max_pauli_weight)
            is_truncated = true
        elseif !isnothing(customtruncfunc) && customtruncfunc(hstr, coeff)
            is_truncated = true
        end

        return is_truncated
    end

    truncate!(truncfunc, prop_cache; kwargs...)
    return 
end 



function PauliPropagation.PropagationBase.truncate!(
    hsum::AbstractHybridSum;
    max_weight::Real=Inf, min_abs_coeff=1e-10,
    max_pauli_weight::Real=Inf, max_majorana_weight::Real=Inf, max_unpaired::Real=Inf,
    max_freq::Real=Inf, max_sins::Real=Inf,
    unpaired_mask=nothing,
    customtruncfunc=nothing,
    fermion_filter=nothing, qubit_filter=nothing, 
    kwargs...
)
    if isnothing(unpaired_mask)
        unpaired_mask = create_unpaired_mask(nfermions(hsum))
    end
    if isnothing(fermion_filter) || isnothing(qubit_filter)
        fermion_filter, qubit_filter = create_filters(hsum)
    end 

    function truncfunc(hstr, coeff)
        is_truncated = false
        if PauliPropagation.truncatemincoeff(coeff, min_abs_coeff)
            is_truncated = true
        elseif MajoranaPropagation.truncateunpaired(typeof(unpaired_mask)(hstr&fermion_filter), max_unpaired, unpaired_mask)
            is_truncated = true
        elseif MajoranaPropagation.truncatemajoranaweight((hstr&fermion_filter), max_majoranan_weight)
            is_truncated = true
        elseif PauliPropagation.truncateweight((hstr&qubit_filter), max_pauli_weight)
            is_truncated = true
        elseif !isnothing(customtruncfunc) && customtruncfunc(hstr, coeff)
            is_truncated = true
        end

        return is_truncated
    end
    hsum = truncate!(truncfunc, hsum; kwargs...)

    return hsum
end