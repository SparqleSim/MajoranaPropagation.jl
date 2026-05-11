
abstract type AbstractHybridSum <: PauliPropagation.AbstractTermSum end

struct HybridSum{TT<:Integer,CT} <: AbstractHybridSum
    nfermionic_sites::Int
    nqubits::Int
    is_spinful::Bool
    terms::Dict{TT,CT}
end

#--- Standard Constructors-------------------------------------
#Initialisation with defined coefficient type
function HybridSum(nfermionic_sites::Int, nqubits::Int, is_spinful::Bool, ::Type{CT}) where {CT}
    if is_spinful
        nfermions = 2 * nfermionic_sites
    else
        nfermions = nfermionic_sites
    end
    TT = getinttype(nqubits + nfermions)
    return HybridSum(nfermionic_sites, nqubits, is_spinful, Dict{TT,CT}())
end

#Default initialisation for empty HybridSumm :
function HybridSum(nfermionic_sites::Int, nqubits::Int, is_spinful::Bool)
    return HybridSum(nfermionic_sites, nqubits, is_spinful, Float64)
end

function HybridSum(msum::MajoranaPropagation.AbstractMajoranaSum, pstr::PauliString)
    nf = MajoranaPropagation.nfermions(msum)
    TT = PauliPropagation.getinttype(nf + pstr.nqubits)
    hybrid_terms = Dict{TT,MajoranaPropagation.coefftype(msum)}()
    hybrid_pstr = TT(pstr.term) << (2 * nf)

    if length(msum.Majoranas) > 0
        for (mstr, coeff) in msum.Majoranas
            hybrid_terms[(TT(mstr)|hybrid_pstr)] = coeff
        end
    else
        hybrid_terms[hybrid_pstr] = 1.0
    end

    return HybridSum(msum.nsites, pstr.nqubits, msum.is_spinful, hybrid_terms)
end

function HybridSum(msum::MajoranaPropagation.AbstractMajoranaSum, psum::PauliPropagation.AbstractPauliSum)
    nfermions = MajoranaPropagation.nfermions(msum)
    TT = PauliPropagation.getinttype(nfermions + psum.nqubits)
    hybrid_terms = Dict{TT,MajoranaPropagation.coefftype(msum)}()

    for (pstr, pcoeff) in psum.terms
        hybrid_pstr = TT(pstr) << (2 * nfermions)
        for (mstr, mcoeff) in msum.Majoranas
            hybrid_terms[(TT(mstr)|hybrid_pstr)] = pcoeff * mcoeff
        end
    end
    return HybridSum(msum.nsites, psum.nqubits, msum.is_spinful, hybrid_terms)
end

function HybridSum(msum::MajoranaPropagation.AbstractMajoranaSum, nqubits::Int, pstr_term::TT) where {TT<:Integer}
    nf = MajoranaPropagation.nfermions(msum)
    dtype = PauliPropagation.getinttype(nf + nqubits)
    hybrid_terms = Dict{dtype,MajoranaPropagation.coefftype(msum)}()
    hybrid_pstr = dtype(pstr_term) << (2 * nf)

    if length(msum.Majoranas) > 0
        for (mstr, coeff) in msum.Majoranas
            hybrid_terms[(dtype(mstr)|hybrid_pstr)] = coeff
        end
    else
        hybrid_terms[hybrid_pstr] = 1.0
    end
    return HybridSum(msum.nsites, nqubits, msum.is_spinful, hybrid_terms)
end

function HybridSum(msum::MajoranaPropagation.AbstractMajoranaSum)
    return HybridSum(msum.nsites, 0, msum.is_spinful, deepcopy(msum.Majoranas))
end

function HybridSum(psum::PauliPropagation.AbstractPauliSum)
    return HybridSum(0, psum.nqubits, false, deepcopy(psum.terms))
end

function HybridSum(pstr::PauliString)
    hybrid_terms = Dict{typeof(pstr.term),typeof(pstr.coeff)}()
    hybrid_terms[pstr.term] = 1.0
    return HybridSum(0, pstr.nqubits, false, hybrid_terms)
end


#--- Constructors for Majorana rotations----------------------------
function HybridSum(n_fermions::Integer, f_symb::Symbol, f_sites::Vector{Int},
    n_qubits::Integer, q_symbols::Vector{Symbol}, q_indices::Vector{Int};
    del_m_id::Bool=false)
    """
    Constructor to be used to identify Majorana rotations per fermionic gate
    """
    TT = getinttype(n_qubits + n_fermions)

    #Single PauliString denoting the Qubit part
    pstr = PauliString(n_qubits + n_fermions, q_symbols, q_indices .+ n_fermions).term

    if isempty(f_sites)
        return HybridSum(n_fermions, n_qubits, false, Dict(pstr => 1.0))
    else
        #Convert the Fermionic gate into MajoranaRotation components
        fermionic_part = MajoranaSum(n_fermions, Val(f_symb), f_sites)

        if del_m_id #remove coefficient associated to identity
            MajoranaPropagation.pop_id!(fermionic_part)
        end

        #Combination of fermionic and qubit parts
        hsum_dict = Dict{TT,MajoranaPropagation.coefftype(fermionic_part)}()
        for (ms, coeff) in fermionic_part
            hs = pstr | convert(TT, ms)
            hsum_dict[hs] = coeff
        end
        return HybridSum(n_fermions, n_qubits, fermionic_part.is_spinful, hsum_dict)
    end
end

#--- Base Overloads ------------------------------------
Base.iterate(hsum::HybridSum, state=1) = iterate(hsum.terms, state)

function Base.length(hsum::HybridSum)
    return length(hsum.terms)
end

function Base.show(io::IO, hsum::HybridSum)
    max_display = 8
    print(io, "HybridSum with $(length(hsum)) term$(length(hsum) == 1 ? "" : "s"):(")
    for (i, (term, coeff)) in enumerate(hsum.terms)
        if i <= max_display
            majorana_string = reverse(string(term; base=2, pad=(2 * hsum.nfermionic_sites)))[1:(2*hsum.nfermionic_sites)]
            pauli_string = ""
            if hsum.nqubits > 0
                pauli_string = inttostring((term >> (2 * hsum.nfermionic_sites)), hsum.nqubits)
            end
            print(io, "\n")
            print(io, "    $(coeff) * $(majorana_string) x $(pauli_string)")
        else
            print(io, "\n    ...")
            break
        end
    end
    print(io, ") \n")
end

function Base.:(==)(hs1::HybridSum, hs2::HybridSum)
    if hs1.nqubits != hs2.nqubits || hs1.nfermionic_sites != hs2.nfermionic_sites || hs1.is_spinful != hs2.is_spinful
        return false
    end
    return hs1.terms == hs2.terms
end

function Base.mergewith!(hsum1::HybridSum{TT,CT}, hsum2::HybridSum{TT,CT}) where {TT<:Integer,CT}
    mergewith!(hsum1.terms, hsum2.terms)
    return hsum1
end

function Base.empty!(hsum::HybridSum{TT,CT}) where {TT<:Integer,CT}
    empty!(hsum.terms)
    return hsum
end

function Base.delete!(hsum::HybridSum{TT,CT}, hs::HybridString{TT}) where {TT<:Integer,CT}
    println("Potentially deprecated method called: delete!(hsum::HybridSum, hs::HybridString)")
    delete!(hsum.terms, hs.term)
end
function Base.delete!(hsum::HybridSum{TT,CT}, hs::TT) where {TT<:Integer,CT}
    delete!(hsum.terms, hs)
end

#--- PropagationBase Overloads -----------------------------------------------
PauliPropagation.PropagationBase.storage(hsum::HybridSum) = hsum.terms
PauliPropagation.PropagationBase.nsites(hsum::HybridSum) = hsum.nfermionic_sites

#--- Other Functions ------------------------------------------------------
is_spinful(hsum::HybridSum) = hsum.is_spinful
nqubits(hsum::HybridSum) = hsum.nqubits

function nfermions(hsum::HybridSum)
    if hsum.is_spinful
        return 2 * hsum.nfermionic_sites
    else
        return hsum.nfermionic_sites
    end
end

function pop_id!(hsum::HybridSum{TT,CT}) where {TT<:Integer,CT}
    if haskey(hsum.terms, 0)
        delete!(hsum.terms, 0)
    end
    return
end

function coefftype(hsum::HybridSum)
    return valtype(collect(values(hsum.terms)))
end

function similar(hsum::HybridSum)
    return HybridSum(hsum.nfermionic_sites, hsum.nqubits, hsum.is_spinful, coefftype(hsum))
end

function set!(hsum::HybridSum{TT,CT}, hs::HybridString{TT}, value::CT) where {TT<:Integer,CT}
    println("Potentially deprecated method called: set!(hsum::HybridSum{TT, CT}, hs::HybridString{TT}, value::CT)")
    set!(hsum, hs.term, value)
    return
end

function set!(hsum::HybridSum{TT,CT}, hs::TT, value::CT) where {TT<:Integer,CT}
    hsum.terms[hs] = value
    return
end

function add!(hsum::HybridSum{TT,CT}, hstr::HybridString{TT}, value::CT) where {TT<:Integer,CT}
    if haskey(hsum.terms, hstr.term)
        hsum.terms[hstr.term] += value
    else
        hsum.terms[hstr.term] = value
    end
    return hsum
end

function add!(hsum::HybridSum{TT,CT}, hs::TT, value::CT) where {TT<:Integer,CT}
    if haskey(hsum.terms, hs)
        hsum.terms[hs] += value
    else
        hsum.terms[hs] = value
    end
    return hsum
end

function mergeandempty!(hsum1::HS, hsum2::HS) where {HS<:AbstractHybridSum}
    mergewith!(hsum1, hsum2)
    empty!(hsum2)
    return hsum1, hsum2
end

#--- Other Stuff -------------------------------------
function create_filters(nfermionic_sites::Int, nqubits::Int, is_spinful::Bool)
    nfermions = is_spinful ? (2 * nfermionic_sites) : nfermionic_sites
    TT = getinttype(nqubits + nfermions)

    #Filter mask for Fermions: 00...0011...11
    fermions_filter = TT(2)^(2 * nfermions) - TT(1)

    #Filter mask for Qubits: 11...1100...00
    qubits_filter = (TT(2)^(2 * nqubits) - TT(1)) << (2 * nfermions)

    #Check edge cases
    if nfermions == 0
        qubits_filter = (TT(2)^(2 * nqubits - 1) - TT(1)) << 1
        qubits_filter |= TT(1)
        fermions_filter = TT(0)
    end
    if nqubits == 0
        fermions_filter = (TT(2)^(2 * nfermions - 1) - TT(1)) << 1
        fermions_filter |= TT(1)
        qubits_filter = TT(0)
    end

    return fermions_filter, qubits_filter
end


function create_filters(hsum::AbstractHybridSum)
    return create_filters(hsum.nfermionic_sites, hsum.nqubits, hsum.is_spinful)
end