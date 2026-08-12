using MajoranaPropagation
using PauliPropagation

"""
    JordanWigner(msum::MajoranaSum)

Map a `MajoranaSum` to a `PauliPropagation.PauliSum` on ``n`` qubits (``n`` the number of fermionic modes) via the Jordan-Wigner transformation ``\\gamma_{2j-1} = Z_1 \\cdots Z_{j-1} X_j``, ``\\gamma_{2j} = Z_1 \\cdots Z_{j-1} Y_j``.
Phases arising from operator reordering are absorbed into the Pauli coefficients.
"""
function JordanWigner(msum::MajoranaSum{TT,CT}) where {TT<:Integer,CT}
    n_fermions = MajoranaPropagation.nfermions(msum)

    psum = PauliSum(CT, n_fermions)
    for (ms, coeff) in msum
        mstr_jw, phase = JordanWigner(ms, msum.nsites, msum.is_spinful)
        add!(psum, mstr_jw, coeff * phase)
    end
    return psum
end

"""
    JordanWigner(mstr::Integer, n_sites::Int, is_spinful::Bool)

Map a single integer-encoded Majorana string to its Jordan-Wigner image.
Returns a tuple `(pstr, phase)` of the integer-encoded Pauli string and the complex phase accumulated by the mapping.
"""
function JordanWigner(mstr::TT, n_sites::Int, is_spinful::Bool) where {TT<:Integer}
    n_fermions = is_spinful ? 2 * n_sites : n_sites
    mstr_jw = TT(0)
    phase = (1im)^omega_L_mult(mstr)
    for i = 1:n_fermions
        gammas = TT(getpauli(mstr, i))
        if gammas == 1 || gammas == 2
            if i > 1
                jw_string = PauliString(n_fermions, repeat([:Z], i - 1), 1:i-1).term
            else
                jw_string = TT(0)
            end
            gammas = jw_string | (gammas << (2 * (i - 1)))
        elseif gammas == 3
            gammas = gammas << (2 * (i - 1))
            phase *= 1im
        end
        mstr_jw, new_sign = pauliprod(mstr_jw, gammas)
        phase *= new_sign
    end
    return mstr_jw, phase
end

function _make_gate(Pj, nq)
    Pj_paulis = inttosymbol(Pj, nq)
    symbs = []
    indices = []
    for (k, p) in enumerate(Pj_paulis)
        if p != :I
            push!(symbs, p)
            push!(indices, k)
        end
    end
    return PauliRotation(symbs, indices)
end

"""
    JordanWigner(n_sites, is_spinful, circ::Vector{FermionicRotation}, thetas::Vector)

Map a circuit of `FermionicRotation` gates with angles `thetas` to an equivalent circuit of `PauliPropagation.PauliRotation` gates.
Returns a tuple `(pp_circ, pp_thetas)`; each fermionic gate expands into one `PauliRotation` per Majorana rotation it contains, with angles rescaled to match `PauliRotation`'s ``e^{-i \\theta/2 P}`` convention.
"""
function JordanWigner(n_sites, is_spinful, circ::Vector{FermionicRotation}, thetas::Vector{CT}) where {CT}
    pp_circ = PauliRotation[]
    pp_thetas = CT[]
    n_fermions = is_spinful ? 2 * n_sites : n_sites
    for (gate, theta) in zip(circ, thetas)
        ms_rotations, coeffs, _ = MajoranaPropagation.getmajoranarotations(gate, n_sites)
        for (ms_rotation, coeff) in zip(ms_rotations, coeffs)
            ms_jw, phase = JordanWigner(ms_rotation.ms_int, n_sites, is_spinful)
            push!(pp_circ, _make_gate(ms_jw, n_fermions))
            # factor 2 because `FermionicRotation` applies each `MajoranaRotation` with angle 2 * coeff * theta,
            # matching `PauliRotation`'s exp(-i * theta/2 * pstr) convention
            push!(pp_thetas, 2.0 * coeff * theta * phase)
        end
    end
    return pp_circ, pp_thetas
end

"""
    majoranapropagation2yao(n_sites, is_spinful, circ, thetas)
    majoranapropagation2yao(msum::MajoranaSum)

Convert a fermionic circuit with angles `thetas`, or a `MajoranaSum`, to a Yao.jl object by first applying the Jordan-Wigner transformation and then `PauliPropagation.paulipropagation2yao`.
Requires Yao.jl to be loaded.
"""
function majoranapropagation2yao(n_sites, is_spinful, circ, thetas)
    nqubits = is_spinful ? 2 * n_sites : n_sites
    pp_circ, pp_thetas = JordanWigner(n_sites, is_spinful, circ, thetas)
    return paulipropagation2yao(nqubits, pp_circ, pp_thetas)
end

function majoranapropagation2yao(msum::MajoranaSum)
    return paulipropagation2yao(JordanWigner(msum))
end
