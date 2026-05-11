
struct HybridRotation{TT<:Integer}
    term::TT
    coeff::Float64
end

struct HybridGate <: PauliPropagation.ParametrizedGate
    rotations::Vector{HybridRotation}
    nfermionic_sites::Int
    nqubits::Int
    is_spinful::Bool
    truncate_each_hrot::Bool
end

#To do: initialisation for non-PauliRotation gates in qubit-part
function HybridGate(nfermionic_sites::Int, f_symbol::Symbol, f_sites::Vector{Int},
    nqubits::Int, q_symbols::Vector{Symbol}, q_sites::Vector{Int},
    is_spinful::Bool)
    nfermions = is_spinful ? (2 * nfermionic_sites) : nfermionic_sites
    rotations::Vector{HybridRotation} = []
    truncate_after_each_hybrot = false

    if nfermions == 0
        pstr = PauliString(nqubits, q_symbols, q_sites)
        push!(rotations, HybridRotation(pstr.term, 1.0))
    elseif nqubits == 0
        fermionic_part = MajoranaSum(nfermions, Val(f_symbol), f_sites)
        for (ms, coeff) in fermionic_part
            push!(rotations, HybridRotation(ms, coeff))
        end
        truncate_after_each_hybrot = !MajoranaPropagation.flag_non_number_preserving(f_symbol)
    else
        hsum = HybridSum(nfermions, f_symbol, f_sites, nqubits, q_symbols, q_sites, del_m_id=true)
        for (hs, coeff) in hsum
            push!(rotations, HybridRotation(hs, coeff))
        end
        truncate_after_each_hybrot = !MajoranaPropagation.flag_non_number_preserving(f_symbol)
    end
    return HybridGate(rotations, nfermionic_sites, nqubits, is_spinful, truncate_after_each_hybrot)
end

#No input for is_spinful => defaultes to false
HybridGate(nfermionic_sites::Int, f_symbol::Symbol, f_sites::Vector{Int}, nqubits::Int, q_symbols::Vector{Symbol}, q_sites::Vector{Int}) = HybridGate(nfermionic_sites, f_symbol, f_sites, nqubits, q_symbols, q_sites, false)
HybridGate(nfermionic_sites::Int, f_symbol::Symbol, f_sites::Vector{Int}, nqubits::Int, is_spinful::Bool) = HybridGate(nfermionic_sites, f_symbol, f_sites, nqubits, Symbol[], Int[], is_spinful)
HybridGate(nfermionic_sites::Int, f_symbol::Symbol, f_sites::Vector{Int}, nqubits::Int) = HybridGate(nfermionic_sites, f_symbol, f_sites, nqubits, Symbol[], Int[], false)
HybridGate(nfermionic_sites::Int, nqubits::Int, q_symbols::Vector{Symbol}, q_sites::Vector{Int}, is_spinful::Bool) = HybridGate(nfermionic_sites, :n, Int[], nqubits, q_symbols, q_sites, is_spinful)
HybridGate(nfermionic_sites::Int, nqubits::Int, q_symbols::Vector{Symbol}, q_sites::Vector{Int}) = HybridGate(nfermionic_sites, :n, Int[], nqubits, q_symbols, q_sites, false)


function _applycos(coeff, cos_theta)
    return coeff * cos_theta
end
function _applysin(coeff, sin_theta)
    return coeff * sin_theta
end


function PauliPropagation.PropagationBase.applytoall!(gate::HybridRotation, prop_cache::HybridPropagationCache, parameter; qubits_filter, fermions_filter, kwargs...)
    #=
    Overload
    Applies one Hybrid rotation to the entire Hybrid sum 
    =#
    hsum = mainsum(prop_cache)
    aux_hsum = auxsum(prop_cache)

    theta = parameter * gate.coeff * 2.0 # multiply coefficient by 2 since exponential implements exp(-i * theta/2 * hstring)
    cos_val = cos(theta)
    sin_val = sin(theta)

    # separate gate into fermionic and qubits part
    P_gate = gate.term & qubits_filter
    mu_gate = gate.term & fermions_filter

    for (hs, coeff) in hsum
        P_string = hs & qubits_filter
        mu_string = hs & fermions_filter

        if PauliPropagation._bitcommutes(P_gate, P_string) == MajoranaPropagation.commutes(mu_gate, mu_string)
            continue
        else
            coeff1 = _applycos(coeff, cos_val)

            sign, new_ms = ms_mult(mu_gate, mu_string, nfermions(hsum))
            Pk, pk_sign = pauliprod(P_string, P_gate)
            coeff2 = _applysin(coeff, sin_val * real(-1im * sign * pk_sign)) #_applysin(coeff, sin_val * -imag(sign * pk_sign))
            hs2 = Pk | new_ms

            set!(hsum, hs, coeff1)
            set!(aux_hsum, hs2, coeff2)
        end
    end
    return
end


function PauliPropagation.PropagationBase.applymergetruncate!(h_gate::HybridGate, prop_cache::AbstractHybridPropagationCache, theta; truncate_each_hr=nothing, kwargs...)
    truncate_after_each_hybrot = h_gate.truncate_each_hrot
    if !isnothing(truncate_each_hr)
        truncate_after_each_hybrot = truncate_each_hr
    end

    # iterate over individual hybrid rotations and apply them to the hybrid sum
    for single_gate in h_gate.rotations
        applytoall!(single_gate, prop_cache, theta; kwargs...)

        # merge the auxiliary Hybrid sum into the original one and empty the auxiliary one
        merge!(prop_cache; kwargs...)

        if truncate_after_each_hybrot
            truncate!(prop_cache; kwargs...)
        end
    end

    if !truncate_after_each_hybrot
        truncate!(prop_cache; kwargs...)
    end

    return prop_cache
end




function PauliPropagation.apply(gate::HybridGate, hstr, coeff, parameter; kwargs...)
    #=
    applies one Hybrid gate to a single string
    =#
    makeshift_hsum = HybridSum(gate.nfermionic_sites, gate.nqubits, gate.is_spinful)
    add!(makeshift_hsum, hstr, coeff)
    applymergetruncate!(gate, makeshift_hsum, similar(makeshift_hsum), [parameter], 1; kwargs...)

    output_tuples = []
    for (h_str, h_coeff) in makeshift_hsum
        push!(output_tuples, (h_str, h_coeff))
    end
    return Tuple(output_tuples)
end