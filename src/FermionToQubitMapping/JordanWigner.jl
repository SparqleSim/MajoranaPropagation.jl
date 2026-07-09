
using MajoranaPropagation
using PauliPropagation


function JordanWigner(msum::MajoranaSum{TT,CT}) where {TT<:Integer,CT}
    @assert msum.is_spinful == false
    n_fermions = MajoranaPropagation.nfermions(msum)

    psum = PauliSum(CT, n_fermions)
    for (ms, coeff) in msum
        mstr_jw, phase = JordanWigner(ms, msum.nsites, msum.is_spinful)
        #@show bitstring(mstr_jw), phase
        add!(psum, mstr_jw, coeff * phase)
    end
    return psum
end

function JordanWigner(mstr::TT, n_sites::Int, is_spinful::Bool) where {TT<:Integer}
    n_fermions = is_spinful ? 2 * n_sites : n_sites
    mstr_jw = TT(0)
    phase = (1im)^omega_L_mult(mstr)
    #@show bitstring(mstr)
    for i = 1:n_fermions
        #println("------")
        gammas = TT(getpauli(mstr, i))
        #@show i, bitstring(gammas)
        if gammas == 1 || gammas == 2
            if i > 1
                jw_string = PauliString(n_fermions, repeat([:Z], i - 1), 1:i-1).term
            else
                jw_string = TT(0)
            end
            #@show bitstring(jw_string)
            gammas = jw_string | (gammas << (2 * (i - 1)))
            #@show bitstring(gammas)
        elseif gammas == 3
            gammas = gammas << (2 * (i - 1))
            phase *= 1im
        end
        mstr_jw, new_sign = pauliprod(mstr_jw, gammas)
        phase *= new_sign
    end
    return mstr_jw, phase
end

function make_gate(Pj, nq)
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

function JordanWigner(circ::Vector{FermionicRotation}, thetas::Vector{CT}, n_sites, is_spinful) where {CT}
    pp_circ = PauliRotation[]
    pp_thetas = CT[]
    for (gate, theta) in zip(circ, thetas)
        ms_rotations, theta_into_majorana_decomposition, _ = MajoranaPropagation.getmajoranarotations(gate, n_sites, theta)
        for (ms_rotation, theta_coeff) in zip(ms_rotations, theta_into_majorana_decomposition)
            ms_jw, phase = JordanWigner(ms_rotation.ms_int, n_sites, is_spinful)
            push!(pp_circ, make_gate(ms_jw, n_sites))
            # factor 2 because `FermionicRotation` applies each `MajoranaRotation` with angle 2 * coeff * theta,
            # matching `PauliRotation`'s exp(-i * theta/2 * pstr) convention
            push!(pp_thetas, 2. * theta_coeff * phase)
        end
    end
    return pp_circ, pp_thetas
end

