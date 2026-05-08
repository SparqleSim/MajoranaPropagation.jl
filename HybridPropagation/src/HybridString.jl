#=
Convention:     - first 2nf bits designate the ferminonic system
                - last 2nq bits designate the qubit system
=#

struct HybridString{TT<:Integer}
    nfermions::Int
    nqubits::Int
    term::TT
end

#Construct from a Pauli and a Majorana string
function HybridString(mString::MajoranaString, pString::PauliString)
    TT = getinttype(mString.nfermions + pString.nqubits)
    hybrid = TT(0)
    hybrid |= mString.gammas
    hybrid |= (TT(pString.term) << (2 * mString.nfermions))
    return HybridString(mString.nfermions, pString.nqubits, hybrid)
end

function nfermions(hs::HybridString)
    return hs.nfermions
end

function nqubits(hs::HybridString)
    return hs.nqubits
end

function Base.show(io::IO, hstr::HybridString)
    majorana_string = reverse(string(hstr.term; base=2, pad=(2 * hstr.nfermions)))[1:(2*hstr.nfermions)]
    pauli_string = inttostring((hstr.term >> (2 * hstr.nfermions)), hstr.nqubits)
    if length(pauli_string) > 20
        pauli_string = pauli_string[1:20] * "..."
    end
    if length(majorana_string) > 20
        majorana_string = majorana_string[1:20] * "..."
    end
    print(io, "HybridString(nfermions: $(hstr.nfermions), nqubits: $(hstr.nqubits), $(majorana_string) x $(pauli_string))")
end

#=
function Base.print(hs::HybridString)
    majorana_string = reverse(string(hstr.term; base=2, pad=(2 * hstr.nfermions)))[1:(2*hstr.nfermions)]
    pauli_string = inttostring((hstr.term >> (2 * hstr.nfermions)), hstr.nqubits)
    print("$(majorana_string) x $(pauli_string)")
end

function Base.println(hs::HybridString)
    print(hs)
    println()
end
=#