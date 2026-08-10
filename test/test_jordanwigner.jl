# Cross-check of the FermionToQubitMappings (Jordan-Wigner) module against the
# independent dense oracle in testutils_dense.jl. The module routes through src
# conventions (omega_L_mult, getmajoranarotations), so agreement with the
# from-scratch dense matrices validates the mapping without circularity.
using MajoranaPropagation
using MajoranaPropagation.FermionToQubitMappings
using PauliPropagation
using Yao
using LinearAlgebra
using Test
using Random
Random.seed!(42)

if !isdefined(Main, :dense_majoranas)
    include("testutils_dense.jl")
end

const _PAULI = Dict(:I => _I2, :X => _X, :Y => _Y, :Z => _Z)

# Dense matrix of a PauliSum in the oracle's kron convention (qubit 1 most
# significant; the JW map sends Majorana mode / qubit index i to qubit i).
function dense_paulisum(psum, nq::Int)
    M = zeros(ComplexF64, 2^nq, 2^nq)
    for (pstr, coeff) in psum
        M .+= coeff .* reduce(kron, [_PAULI[s] for s in inttosymbol(pstr, nq)])
    end
    return M
end

function dense_rotation_generator(g::PauliRotation, nq::Int)
    ops = Matrix{ComplexF64}[_I2 for _ in 1:nq]
    for (s, q) in zip(g.symbols, g.qinds)
        ops[q] = _PAULI[s]
    end
    return reduce(kron, ops)
end

function random_odd_string(nf::Int)
    g = random_string(nf)
    while get_weight(g) % 2 != 1
        g = random_string(nf)
    end
    return g
end

@testset "Jordan-Wigner module vs dense oracle" begin

    @testset "operator mapping, spinless" begin
        nf = 4
        gam = dense_majoranas(nf)
        for complexcoeffs in (false, true), _ in 1:10
            msum = random_even_msum(nf, 5; complexcoeffs)
            # include an odd-weight (unpaired) string: the Z-string handling
            # in the map is only exercised by these
            PauliPropagation.PropagationBase.add!(msum, random_odd_string(nf), randn())
            @test isapprox(dense_paulisum(JordanWigner(msum), nf), dense_sum(msum, gam), atol=1e-12)
        end

        # identity string maps to the identity with the coefficient untouched
        msum = MajoranaSum(Float64, nf)
        PauliPropagation.PropagationBase.add!(msum, getinttype(nf)(0), 0.75)
        @test isapprox(dense_paulisum(JordanWigner(msum), nf), 0.75 * LinearAlgebra.I(2^nf), atol=1e-14)

        # human-readable pins: n_2 -> (I - Z_2)/2 and adjacent hopping
        # (Z-strings cancel) -> (X_2 X_3 + Y_2 Y_3)/2
        @test isapprox(dense_paulisum(JordanWigner(MajoranaSum(nf, :n, 2)), nf),
            (kron(_I2, _I2, _I2, _I2) - kron(_I2, _Z, _I2, _I2)) / 2, atol=1e-14)
        @test isapprox(dense_paulisum(JordanWigner(MajoranaSum(nf, :hop, [2, 3])), nf),
            (kron(_I2, _X, _X, _I2) + kron(_I2, _Y, _Y, _I2)) / 2, atol=1e-14)
    end

    @testset "operator mapping, spinful" begin
        n_sites = 2
        nf = 2 * n_sites
        gam = dense_majoranas(nf)
        for _ in 1:10
            msum = random_even_msum(nf, 5; is_spinful=true)
            @test isapprox(dense_paulisum(JordanWigner(msum), nf), dense_sum(msum, gam), atol=1e-12)
        end
        # site s: up spin on qubit 2s-1, down spin on qubit 2s
        @test isapprox(dense_paulisum(JordanWigner(MajoranaSum(n_sites, :nup, 2)), nf),
            (kron(_I2, _I2, _I2, _I2) - kron(_I2, _I2, _Z, _I2)) / 2, atol=1e-14)
        @test isapprox(dense_paulisum(JordanWigner(MajoranaSum(n_sites, :ndn, 1)), nf),
            (kron(_I2, _Z, _I2, _I2) .* -1 + kron(_I2, _I2, _I2, _I2)) / 2, atol=1e-14)
    end

    @testset "circuit translation" begin
        # the translated Pauli circuit's unitary must equal the product of
        # exp(-i theta op) fermionic unitaries exactly (no Trotter slack:
        # the expansion of each FermionicRotation commutes internally)
        for (n_sites, is_spinful, gatespecs) in (
            (4, false, [(:hop, [1, 2]), (:nn, [2, 4]), (:hop, [3, 4])]),
            (2, true, [(:hopup, [1, 2]), (:hopdn, [1, 2]), (:nupndn, 1)]),
        )
            nf = is_spinful ? 2 * n_sites : n_sites
            gam = dense_majoranas(nf)
            circ = [FermionicRotation(s, q) for (s, q) in gatespecs]
            thetas = randn(length(circ))

            pp_circ, pp_thetas = JordanWigner(n_sites, is_spinful, circ, thetas)

            U_pauli = Matrix{ComplexF64}(LinearAlgebra.I, 2^nf, 2^nf)
            for (g, th) in zip(pp_circ, pp_thetas)
                U_pauli = exp(-im * th / 2 * dense_rotation_generator(g, nf)) * U_pauli
            end
            U_ferm = Matrix{ComplexF64}(LinearAlgebra.I, 2^nf, 2^nf)
            for ((s, q), th) in zip(gatespecs, thetas)
                op = MajoranaSum(n_sites, s, q)
                # the rotation decomposition drops the identity component of
                # e.g. :nn (it only contributes a global phase e^{-i theta c0}),
                # so drop it on the dense side too to compare exactly
                MajoranaPropagation.pop_id!(op)
                U_ferm = exp(-im * th * dense_sum(op, gam)) * U_ferm
            end
            @test isapprox(U_pauli, U_ferm, atol=1e-12)
        end
    end

    @testset "Yao integration" begin
        # expectation values of the mapped Yao operator on every Fock product
        # state agree with the dense oracle
        nf = 3
        gam = dense_majoranas(nf)
        for _ in 1:5
            msum = random_even_msum(nf, 4)
            yao_obs = majoranapropagation2yao(msum)
            dM = dense_sum(msum, gam)
            for f in allfockstates(nf)
                occupied = [s for s in 1:nf if _fock_occupied(f, s)]
                psi = zero_state(nf)
                if !isempty(occupied)
                    Yao.apply!(psi, chain(nf, put(q => Yao.X) for q in occupied))
                end
                v = dense_fock(f)
                @test isapprox(Yao.expect(yao_obs, psi), v' * dM * v, atol=1e-12)
            end
        end
    end
end
