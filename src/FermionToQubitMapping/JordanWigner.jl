
using MajoranaPropagation
using PauliPropagation

function JordanWigner(msum::MajoranaSum{TT,CT}) where {TT<:Integer,CT}
    n_fermions = MajoranaPropagation.nfermions(msum)

    psum = PauliSum(CT, n_fermions)
    for (ms, coeff) in msum
        mstr_jw, phase = JordanWigner(ms, msum.nsites, msum.is_spinful)
        add!(psum, mstr_jw, coeff * phase)
    end
    return psum
end

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

function majoranapropagation2yao(n_sites, is_spinful, circ, thetas)
    nqubits = is_spinful ? 2 * n_sites : n_sites
    pp_circ, pp_thetas = JordanWigner(n_sites, is_spinful, circ, thetas)
    return paulipropagation2yao(nqubits, pp_circ, pp_thetas)
end

function majoranapropagation2yao(msum::MajoranaSum)
    return paulipropagation2yao(JordanWigner(msum))
end
