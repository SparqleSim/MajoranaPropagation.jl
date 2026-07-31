
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
    for (gate, theta) in zip(circ, thetas)
        ms_rotations, coeffs, _ = MajoranaPropagation.getmajoranarotations(gate, n_sites)
        for (ms_rotation, coeff) in zip(ms_rotations, coeffs)
            # TODO: buggy if spinful, FIX asap
            ms_jw, phase = JordanWigner(ms_rotation.ms_int, n_sites, is_spinful)
            push!(pp_circ, _make_gate(ms_jw, n_sites))
            # factor 2 because `FermionicRotation` applies each `MajoranaRotation` with angle 2 * coeff * theta,
            # matching `PauliRotation`'s exp(-i * theta/2 * pstr) convention
            push!(pp_thetas, 2.0 * coeff * theta * phase)
        end
    end
    return pp_circ, pp_thetas
end

# PauliPropagation >= 0.8 defines `paulipropagation2yao` itself; with older versions
# it only exists in YaoBlocks' PauliPropagationExt, so resolve it at call time
function _paulipropagation2yao(args...)
    if isdefined(PauliPropagation, :paulipropagation2yao)
        return PauliPropagation.paulipropagation2yao(args...)
    end
    for (id, mod) in Base.loaded_modules
        if id.name == "YaoBlocks" && isdefined(mod, :paulipropagation2yao)
            return mod.paulipropagation2yao(args...)
        end
    end
    error("`paulipropagation2yao` not found. Load Yao (`using Yao`) or use PauliPropagation >= 0.8.")
end

function majoranapropagation2yao(n_sites, is_spinful, circ, thetas)
    nqubits = is_spinful ? 2 * n_sites : n_sites
    pp_circ, pp_thetas = JordanWigner(n_sites, is_spinful, circ, thetas)
    return _paulipropagation2yao(nqubits, pp_circ, pp_thetas)
end

function majoranapropagation2yao(msum::MajoranaSum)
    return _paulipropagation2yao(JordanWigner(msum))
end
