#--- HybridGate ---------------------------------------------------------------------------
"""
HybridGate (G): FermionicGate X PauliGate
The HybridGate is decomposed into multiple HybridRotation objects during propagate!
(Since the FermionicGate needs to be decomposed into MajoranaRotations)
G = Σ_j c_j * R_j
--> exp(-i * Θ * G) = Π_j exp(-i * c_j * Θ * R_j) = Π_j exp(-i * Θ_j * R_j), with Θ_j = c_j * Θ

"""
#!!!!!!!!!!To Do: Implement for pauli gates other than Pauli Rotation
struct HybridGate <: PauliPropagation.ParametrizedGate 
    nfermionic_sites::Int
    f_symbol::Symbol
    f_sites::Vector{Int}
    nqubits::Int
    q_symbols::Vector{Symbol}
    q_sites::Vector{Int}
    is_spinful::Bool
    truncate_each_hrot::Bool 

    function HybridGate(nfermionic_sites::Int, f_symbol::Symbol, f_sites::Vector{Int},
        nqubits::Int, q_symbols::Vector{Symbol}, q_sites::Vector{Int}, 
        is_spinful::Bool)

        count((x -> (x > nqubits)), q_sites) > 0 ? throw(AssertionError("Qubit indices in list out of bounds")) : nothing
        count((x -> (x > nfermionic_sites)), f_sites) > 0 ? throw(AssertionError("Fermion indices in list out of bounds")) : nothing

        #Make use of preimplemented checks
        q_part = PauliRotation(q_symbols, q_sites)

        #Identify truncation scheme
        truncate_individual_HybRot = false
        if nfermionic_sites > 0
            truncate_individual_HybRot = !MajoranaPropagation.flag_non_number_preserving(f_symbol)
        end 

        return new(nfermionic_sites, f_symbol, f_sites, nqubits, q_part.symbols, q_part.qinds, is_spinful, truncate_individual_HybRot)
    end
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

#--- HybridRotations ---------------------------------------------------------------------------
"""
HybridRotation (R): MajoranaRotation X PauliGate
Entity on which applytoall is actually performed on (checking commutation relations & corresponding branching)
--> exp(-i * c_j * Θ * R_j)
term = R_j  (<:Integer)
coeff = c_j
Branching of HybridString S corresponds to: cos(2 * c_j * Θ) * S + i * sin(2 * c_j * Θ) * R_j * S
(Factor 2 does needn't be accounted for in the definition of Θ)
"""
struct HybridRotation{TT<:Integer}
    term::TT
    coeff::Float64
end



#--- Base Overloads ----------------------------------------------------------------------------------
function Base.show(io::IO, gate::HybridGate)
    f_sym = ""
    if !(isempty(gate.f_sites))
        f_sym = f_sym * string(gate.f_symbol) * "_" * string(gate.f_sites)
    end 
    q_sym = ""
    if !(isempty(gate.q_sites))
        for (i, site) in enumerate(gate.q_sites)
            q_sym = q_sym * string(gate.q_symbols[i]) * "_" * string(site) * " "
        end 
    end 
    output = "HybridGate: " * f_sym * " x " * q_sym * "\n"
    print(io, output)
end 


#--- Functions ------------------------------------------------------------------------------------------
"""
Decomposes the fermionic gate in majorana rotations and fuses them with the pauli part
Outputs an array of HybridRotations
"""
function gethybridrotations(gate::HybridGate)
    nfermions = gate.is_spinful ? (2 * gate.nfermionic_sites) : gate.nfermionic_sites
    rotations::Vector{HybridRotation} = []
    if nfermions == 0
        pstr = symboltoint(gate.nqubits, gate.q_symbols, gate.q_sites)
        push!(rotations, HybridRotation(pstr, 1.0))
    elseif gate.nqubits == 0
        for (mstr, coeff) in MajoranaSum(nfermions, Val(gate.f_symbol), gate.f_sites)
            push!(rotations, HybridRotation(mstr, coeff))
        end 
    else 
        hsum = HybridSum(nfermions, gate.f_symbol, gate.f_sites, gate.nqubits, gate.q_symbols, gate.q_sites, del_m_id=true)
        for (hstr, coeff) in hsum
            push!(rotations, HybridRotation(hstr, coeff))
        end 
    end 
    return rotations
end



#--- PropagationBase Overloads ----------------------------------------------------------------------------------
"""
Applies one HybridRotation to the entire PropagationCache
"""
function PauliPropagation.PropagationBase.applytoall!(h_rot::HybridRotation, prop_cache::HybridPropagationCache, parameter; qubits_filter, fermions_filter, kwargs...)
    
    hsum = mainsum(prop_cache)
    aux_hsum = auxsum(prop_cache)

    theta = parameter * h_rot.coeff * 2.0 # multiply coefficient by 2 since exponential implements exp(-i * theta * HybridGate)
    cos_val = cos(theta)
    sin_val = sin(theta)

    # separate gate into fermionic and qubits part
    P_gate = h_rot.term & qubits_filter
    mu_gate = h_rot.term & fermions_filter

    for (hs, coeff) in hsum
        P_string = hs & qubits_filter
        mu_string = hs & fermions_filter

        if PauliPropagation._bitcommutes(P_gate, P_string) == MajoranaPropagation.commutes(mu_gate, mu_string)
            continue
        else
            coeff1 = _applycos(coeff, cos_val)

            sign, new_ms = ms_mult(mu_gate, mu_string, nfermions(hsum))
            Pk, pk_sign = pauliprod(P_gate, P_string)
            coeff2 = _applysin(coeff, sin_val * -imag(sign * pk_sign))#* real(-1im * sign * pk_sign))
            hs2 = Pk | new_ms

            set!(hsum, hs, coeff1)
            set!(aux_hsum, hs2, coeff2)
        end
    end
    return
end



"""
Applies one HybridGate to one HybridGate
"""
function PauliPropagation.PropagationBase.applymergetruncate!(h_gate::HybridGate, prop_cache::AbstractHybridPropagationCache, theta; truncate_each_hr=nothing, kwargs...)

    truncate_after_each_hybrot = h_gate.truncate_each_hrot
    if !isnothing(truncate_each_hr)
        truncate_after_each_hybrot = truncate_each_hr
    end

    hybrid_rotations = gethybridrotations(h_gate)

    # iterate over individual hybrid rotations and apply them to the hybrid sum
    for rotation in hybrid_rotations
        applytoall!(rotation, prop_cache, theta; kwargs...)

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



"""
Applies one HybridGate to a single HybridString term
"""
function PauliPropagation.apply(gate::HybridGate, hstr, coeff, parameter; kwargs...)
    
    makeshift_hsum = HybridSum(gate.nfermionic_sites, gate.nqubits, gate.is_spinful)
    add!(makeshift_hsum, hstr, coeff)
    applymergetruncate!(gate, makeshift_hsum, similar(makeshift_hsum), [parameter], 1; kwargs...)

    output_tuples = []
    for (h_str, h_coeff) in makeshift_hsum
        push!(output_tuples, (h_str, h_coeff))
    end
    return Tuple(output_tuples)
end



#--- MaskedHybridGate --------------------------------------------------------------
"""
Gate struct that already contains the HybridRotations (they don't need to be recomputed each propagation step)
Should therefore be faster
"""
struct MaskedHybridGate <: PauliPropagation.ParametrizedGate
    rotations::Vector{HybridRotation}
    nfermionic_sites::Int
    nqubits::Int
    is_spinful::Bool
    truncate_each_hrot::Bool
end

function MaskedHybridGate(nfermionic_sites::Int, f_symbol::Symbol, f_sites::Vector{Int},
    nqubits::Int, q_symbols::Vector{Symbol}, q_sites::Vector{Int},
    is_spinful::Bool)

    h_gate = HybridGate(nfermionic_sites, f_symbol, f_sites, nqubits, q_symbols, q_sites, is_spinful)
    return MaskedHybridGate(gethybridrotations(h_gate), nfermionic_sites, nqubits, is_spinful, h_gate.truncate_each_hrot)
end 


MaskedHybridGate(nfermionic_sites::Int, f_symbol::Symbol, f_sites::Vector{Int}, nqubits::Int, q_symbols::Vector{Symbol}, q_sites::Vector{Int}) = MaskedHybridGate(nfermionic_sites, f_symbol, f_sites, nqubits, q_symbols, q_sites, false)
MaskedHybridGate(nfermionic_sites::Int, f_symbol::Symbol, f_sites::Vector{Int}, nqubits::Int, is_spinful::Bool) = MaskedHybridGate(nfermionic_sites, f_symbol, f_sites, nqubits, Symbol[], Int[], is_spinful)
MaskedHybridGate(nfermionic_sites::Int, f_symbol::Symbol, f_sites::Vector{Int}, nqubits::Int) = MaskedHybridGate(nfermionic_sites, f_symbol, f_sites, nqubits, Symbol[], Int[], false)
MaskedHybridGate(nfermionic_sites::Int, nqubits::Int, q_symbols::Vector{Symbol}, q_sites::Vector{Int}, is_spinful::Bool) = MaskedHybridGate(nfermionic_sites, :n, Int[], nqubits, q_symbols, q_sites, is_spinful)
MaskedHybridGate(nfermionic_sites::Int, nqubits::Int, q_symbols::Vector{Symbol}, q_sites::Vector{Int}) = MaskedHybridGate(nfermionic_sites, :n, Int[], nqubits, q_symbols, q_sites, false)


"""
Applies one MaskedHybridGate to one HybridGate
"""
function PauliPropagation.PropagationBase.applymergetruncate!(h_gate::MaskedHybridGate, prop_cache::AbstractHybridPropagationCache, theta; truncate_each_hr=nothing, kwargs...)

    truncate_after_each_hybrot = h_gate.truncate_each_hrot
    if !isnothing(truncate_each_hr)
        truncate_after_each_hybrot = truncate_each_hr
    end

    # iterate over individual hybrid rotations and apply them to the hybrid sum
    for rotation in h_gate.rotations
        applytoall!(rotation, prop_cache, theta; kwargs...)

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