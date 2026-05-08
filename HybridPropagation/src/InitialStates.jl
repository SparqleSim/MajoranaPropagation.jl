#=
Basis for the qubit-only EigenState:
- x-basis:  Z term (11) is placed at the -1 eigenstate sites
            sorter is _countbityz
            identifier is _countbity
- y-basis:  X term (01) is placed at the -1 eigenstate sites
            sorter is _countbitxz
            identifier is _countbitz
- z-basis:  Y term (10) is placed at the -1 eigenstate sites
            sorter is _countbityz
            identifier is _countbitx

Terminology:    - sorter: Function that is responsible to "sort out" (return 0 for) the given string in the sum, if the specified
                    Pauli weight its > 0
                - identifier: Function that is responsible to identify the number in the overlapps in -1 Eigenstate sites 
                    and Eigenvalue operations
=#

const basis_symbols = (:x, :y, :z)

function evaluate_basis(basis::Symbol)
    return findfirst(s -> s == basis, basis_symbols)
end

function assign_transform_pauli(basis::Symbol)
    return  1 + mod((evaluate_basis(basis) - 2), 3)
end 

struct EigenState{TT<:Integer}
    n_qubits::Int
    basis::Symbol
    occupied_sites::TT
end 

function EigenState(nqubits::Int, site_list, basis::Symbol; inverted_list::Bool=false)
    (site_list isa Array) || (site_list isa Tuple) ? nothing : throw(AssertionError("Site list needs to be tuple or array"))
    basis in basis_symbols ? nothing : throw(AssertionError("Invalid basis input"))

    TT = PauliPropagation.getinttype(nqubits)
    
    if inverted_list #Construct the mirrored list
        new_list = Int[]
        for j in 1:nqubits
            if !(j in site_list)
                push!(new_list, j)
            end 
        end
        site_list = new_list
    end 

    state = TT(0)
    pauli = TT(assign_transform_pauli(basis))
    for site in site_list
        state |= pauli << (2 * (site - 1))
    end

    return EigenState(nqubits, basis, state)
end

function assign_sortout_func(basis::Symbol)
    if basis == :x
        return PauliPropagation._countbityz
    elseif basis == :y 
        return PauliPropagation._countbitxz
    elseif basis == :z 
        return PauliPropagation._countbitxy 
    else
        throw(AssertionError("Invalid basis input"))
    end
end 

function assign_identifier_func(basis::Symbol)
    if basis == :x
        return PauliPropagation._countbity
    elseif basis == :y 
        return PauliPropagation._countbitz
    elseif basis == :z 
        return PauliPropagation._countbitx
    else
        throw(AssertionError("Invalid basis input"))
    end
end

function overlapwitheigenstate(pstr::TT, state::EigenState, sorter::Function, identifier::Function) where {TT<:Integer}
   """
   pstr contains only the Qubit part
   """ 
    if sorter(pstr) > 0
        return 0
    else
        return (-1)^(identifier(pstr ⊻ state.occupied_sites))
    end 
end 

function Base.show(io::IO, state::EigenState)
    qubit_part = inttostring((state.occupied_sites), state.n_qubits)
    basis = "x"
    if state.basis == :y 
        basis ="y"
    elseif state.basis == :z 
        basis = "z"
    end 

    print(io, "Eigen State, $(state.n_qubits) qubits in ", basis, "-basis: ", qubit_part, "\n")
end 



#--- Hybrid EigenState-------------------------------------------------------------
struct HybridEigenState
    f_part::FockState
    q_part::EigenState

    function HybridEigenState(fermionic_part::FockState, qubit_part::EigenState)
        TT = PauliPropagation.getinttype(fermionic_part.n_sites + qubit_part.n_qubits)
        return new(fermionic_part,  
                    EigenState(qubit_part.n_qubits, qubit_part.basis, (TT(qubit_part.occupied_sites) << (2 * fermionic_part.n_sites))))
    end 
end

function overlapwithstate(observable::AbstractHybridSum, state::HybridEigenState, 
                            fermion_filter::TT, qubit_filter::TT, 
                            unpaired_mask::AT, 
                            sorter::Function, identifier::Function) where {TT<:Integer, AT<:Integer}
    res = 0. 
    for (hstr, coeff) in observable
        phase = overlapwithfock(typeof(unpaired_mask)(fermion_filter & hstr), unpaired_mask, state.f_part)
        if phase != 0 
            phase *= overlapwitheigenstate((qubit_filter & hstr), state.q_part, sorter, identifier)
        end
        res += tonumber(coeff) * phase 
    end 
    return res
end 

function overlapwithstate(observable::AbstractHybridSum, state::HybridEigenState)
    sorter = (pstr) -> assign_sortout_func(state.q_part.basis)(pstr)
    identifier = (pstr) -> assign_identifier_func(state.q_part.basis)(pstr)

    f_filter, q_filter = create_filters(observable)

    unpaired_mask = create_unpaired_mask(nfermions(observable))
    
    return overlapwithstate(observable, state, f_filter, q_filter, unpaired_mask, sorter, identifier)
end 


function Base.show(io::IO, state::HybridEigenState)
    qubit_part = inttostring((state.q_part.occupied_sites >> (2 * state.f_part.n_sites)), state.q_part.n_qubits)

    is_spinful = state.f_part.is_spinful
    fermion_part = "( "
    nf = state.f_part.n_sites

    if is_spinful
        nf *= 2
        for site in 1:(state.f_part.n_sites - 1 )
            fermion_part = fermion_part * string((state.f_part.occupied_sites >> (4 * site - 4)) & 1) * "  " * string((state.f_part.occupied_sites >> (4 * site - 2)) & 1) * " , "
        end 
        if state.f_part.n_sites > 0
            fermion_part = fermion_part * string((state.f_part.occupied_sites >> (4 * state.f_part.n_sites - 4)) & 1) * "  " * string((state.f_part.occupied_sites >> (4 * state.f_part.n_sites - 2)) & 1)
        end 
    else
        for site in 1:(state.f_part.n_sites - 1 )
            fermion_part = fermion_part * string((state.f_part.occupied_sites >> (2 * site - 2)) & 1) * " , "
        end 
        if state.f_part.n_sites > 0
            fermion_part = fermion_part * string((state.f_part.occupied_sites >> (2 * state.f_part.n_sites - 2)) & 1)
        end 
    end 
    fermion_part = fermion_part * " )"
    
    print(io, "Hybrid State: $(nf) fermions ($(state.f_part.is_spinful ? "spinfull" : "spinless")), $(state.q_part.n_qubits) qubits ($(state.q_part.basis)-basis) \n   ", fermion_part, " x ", qubit_part, "\n")
    return 
end