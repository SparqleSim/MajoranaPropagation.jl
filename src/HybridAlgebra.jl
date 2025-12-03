using LinearAlgebra
using Bits

struct hybridString{TT<:Integer}
    nqubits::Int
    nfermions::Int
    string::TT
end

#=function hybridString(nqubits::Int, nfermions::Int, q_symbs, q_inds, f_gammas)
    string = PauliString(nqubits + nfermions, q_symbs, q_inds).term
    TT = typeof(string)
    for i in f_gammas
        string = string | (TT(1) << (i - 1 + 2 * nqubits))
    end
    return hybridString(nqubits, nfermions, string)
end=#

function Base.print(hs::hybridString)
    print("$(reverse(string(hs.string; base=2, pad=2 * (hs.nqubits + hs.nfermions))))")
end
function Base.println(hs::hybridString)
    print(hs)
    println()
end

struct hybridSum{TT<:Integer, CT}
    nqubits::Int
    nfermionic_sites::Int
    is_spinful::Bool
    strings::Dict{TT, CT}
end

function hybridSum(nqubits::Int, nfermionic_sites::Int, is_spinful::Bool, ::Type{CT}) where {CT}
    if is_spinful 
        nfermions = 2 * nfermionic_sites
    else 
        nfermions = nfermionic_sites
    end
    TT = getinttype(nqubits + nfermions)
    return hybridSum(nqubits, nfermionic_sites, is_spinful, Dict{TT,CT}())
end

function nfermions(hsum::hybridSum)
    if hsum.is_spinful
        return 2 * hsum.nfermionic_sites
    end
    return hsum.nfermionic_sites
end

Base.iterate(hsum::hybridSum, state=1) = iterate(hsum.strings, state)

function Base.length(hsum::hybridSum)
    return length(hsum.strings)
end

function Base.show(io::IO, hsum::hybridSum)
    max_display = 8
    print(io, "hybridSum with $(length(hsum)) term$(length(hsum) == 1 ? "" : "s"):(")
    for (i, (hs_value, coeff)) in enumerate(hsum.strings)
        if i <= max_display
            print(io, "\n")
            print(io, "    $(coeff) * $(reverse(string(hs_value; base=2, pad=2 * (hsum.nqubits + nfermions(hsum)))))")
        else
            print(io, "\n    ...")
            break
        end
    end
    print(io, ")")
end

function Base.:(==)(hs1::hybridSum, hs2::hybridSum)
    if hs1.nqubits != hs2.nqubits || hs1.nfermionic_sites != hs2.nfermionic_sites || hs1.is_spinful != hs2.is_spinful
        return false
    end
    return hs1.strings == hs2.strings
end

function coefftype(hsum::hybridSum)
    return valtype(collect(values(hsum.strings)))
end

function similar(hsum::hybridSum)
    return hybridSum(hsum.nqubits, hsum.nfermionic_sites, hsum.is_spinful, coefftype(hsum))
end

function set!(hsum::hybridSum{TT, CT}, hs::hybridString{TT}, value::CT) where {TT<:Integer, CT}
    println("Potentially deprecated method called: set!(hsum::hybridSum{TT, CT}, hs::hybridString{TT}, value::CT)")
    set!(hsum, hs.string, value)
    return
end

function set!(hsum::hybridSum{TT, CT}, hs::TT, value::CT) where {TT<:Integer, CT}
    hsum.strings[hs] = value
    return
end

function sum_add!(hsum::hybridSum{TT, CT}, hs::hybridString{TT}, value::CT) where {TT<:Integer, CT}
    println("Potentially deprecated method called: sum_add!(hsum::hybridSum{TT, CT}, hs::hybridString{TT}, value::CT)")
    sum_add!(hsum, hs.string, value)
end

function sum_add!(hsum::hybridSum{TT, CT}, hs::TT, value::CT) where {TT<:Integer, CT}
    if haskey(hsum.strings, hs)
        hsum.strings[hs] += value 
    else 
        hsum.strings[hs] = value 
    end
end

function Base.mergewith!(merge, hsum1::hybridSum{TT, CT}, hsum2::hybridSum{TT, CT}) where {TT<:Integer, CT}
    mergewith!(merge, hsum1.strings, hsum2.strings)
    return hsum1
end

function Base.empty!(hsum::hybridSum{TT, CT}) where {TT<:Integer, CT}
    empty!(hsum.strings)
    return hsum
end

function Base.delete!(hsum::hybridSum{TT, CT}, hs::hybridString{TT}) where {TT<:Integer, CT}
    println("Potentially deprecated method called: delete!(hsum::hybridSum, hs::hybridString)")
    delete!(hsum.strings, hs.string)
end
function Base.delete!(hsum::hybridSum{TT, CT}, hs::TT) where {TT<:Integer, CT}
    delete!(hsum.strings, hs)
end

function pop_id!(msum::hybridSum{TT, CT}) where {TT<:Integer, CT}
    if haskey(msum.strings, 0)
        delete!(msum.strings, 0)
    end
    return
end

function create_filters(nqubits::Int, n_fermionic_sites::Int, is_spinful::Bool)
    nfermions = is_spinful ? 2 * n_fermionic_sites : n_fermionic_sites
    TT = getinttype(nqubits + nfermions)

    qubits_filter = TT(0)
    for i=1:2*nqubits
        qubits_filter |= TT(1) << (i - 1 + 2 * nfermions)
    end
    println("$(reverse(string(qubits_filter; base=2, pad=2 * (nqubits + nfermions))))")

    fermions_filter = TT(0)
    for i=1:2*nfermions
        fermions_filter |= TT(1) << (i - 1)
    end
    println("$(reverse(string(fermions_filter; base=2, pad=2 * (nqubits + nfermions))))")

    fermions_max_single = create_max_single_filter(nfermions+nqubits) & fermions_filter
    #println("$(reverse(string(fermions_max_single; base=2, pad=2 * (nqubits + nfermions))))")

    return qubits_filter, fermions_filter, fermions_max_single
end

function hs_mult(hsum1::hybridSum{TT, CT}, hsum2::hybridSum{TT, CT}, qubits_filter::TT, fermions_filter::TT) where {TT<:Integer, CT}
    res = similar(hsum1)
    for (hs1, coeff1) in hsum1.strings
        P1::TT = hs1 & qubits_filter
        mu1::TT = hs1 & fermions_filter
        for (hs2, coeff2) in hsum2.strings
            P2::TT = hs2 & qubits_filter
            mu2::TT = hs2 & fermions_filter

            mu_sign, new_ms = ms_mult(mu1, mu2, nfermions(hsum1))
            new_P, P_sign = pauliprod(P1, P2)

            hs_res = new_ms | new_P
            sum_add!(res, hs_res, coeff1 * coeff2 * mu_sign * P_sign)
        end
    end
    return res
end

function Base.:(*)(hsum1::hybridSum{TT, CT}, hsum2::hybridSum{TT, CT}) where {TT<:Integer, CT}
    qubits_filter, fermions_filter, _ = create_filters(hsum1.nqubits, nfermions(hsum1), hsum1.is_spinful)
    return hs_mult(hsum1, hsum2, qubits_filter, fermions_filter)
end

function fock_computational_evaluate(hsum::hybridSum{TT, CT}, fermionic_sites_with_particle, qubit_computational) where {TT<:Integer, CT}
    #TODO: fix 
    res = 0.

    #TT = getinttype(hsum.nqubits + nfermions(hsum))
    TT_q = getinttype(hsum.nqubits)
    TT_f = getinttype(hsum.nfermions)

    n_qubit_bits = 2 * hsum.nqubits
    qubit_mask = (one(TT) << n_qubit_bits) - one(TT)
    get_qubit_part = str -> TT_q(str & qubit_mask)
    
    n_fermion_bits = 2 * hsum.nfermions
    get_fermionic_part = str -> TT_f((str >> n_qubit_bits) & ((one(typeof(str)) << n_fermion_bits) - one(typeof(str))))
    singles_filter = create_max_single_filter(hsum.nfermions)

    for (hs, coeff) in hsum 
        Pi = PauliString(hsum.nqubits, get_qubit_part(hs), 1.)
        mu_v = get_fermionic_part(hs)

        res += coeff * overlapwithcomputational(Pi, qubit_computational) * fockevaluate(mu_v, 1., singles_filter, fermionic_sites_with_particle)
    end
    return res 
end