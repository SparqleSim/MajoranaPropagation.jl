# Fock-state overlaps (diagonal, off-diagonal, superpositions) against the
# dense Jordan-Wigner oracle. This is the CI regression suite for the
# superposition/matrix-element code fixed in issue #15, replacing the
# ITensor-based test_fock.jl which is not part of the test suite.
using MajoranaPropagation
using PauliPropagation
using LinearAlgebra
using Test
using Random
Random.seed!(42)

if !isdefined(Main, :dense_majoranas)
    include("testutils_dense.jl")
end

@testset "Fock overlaps vs dense oracle" begin

    @testset "oracle self-check: gamma algebra" begin
        nf = 4
        gam = dense_majoranas(nf)
        # {gamma_i, gamma_j} = 2 delta_ij
        for i in 1:2nf, j in 1:2nf
            @test isapprox(gam[i] * gam[j] + gam[j] * gam[i], (i == j ? 2.0 : 0.0) * I, atol=1e-13)
        end
        # the number operator from the oracle's a's matches the stored :n convention
        as = dense_annihilators(nf, gam)
        for s in 1:nf
            @test isapprox(dense_sum(MajoranaSum(nf, :n, s), gam), as[s]' * as[s], atol=1e-13)
        end
    end

    @testset "spinless diagonal <f|msum|f>" begin
        nf = 4
        gam = dense_majoranas(nf)
        focks = allfockstates(nf)
        for _ in 1:25
            msum = random_even_msum(nf, 6)
            D = dense_sum(msum, gam)
            for f in focks
                v = dense_fock(f)
                @test abs(overlapwithfock(msum, f) - real(v' * D * v)) < 1e-12
            end
        end

        # dict- and vector-backed sums agree
        msum = random_even_msum(nf, 6)
        vsum = VectorMajoranaSum(msum)
        for f in focks
            @test overlapwithfock(msum, f) == overlapwithfock(vsum, f)
        end
    end

    @testset "diagonal edge cases" begin
        nf = 3
        f = FockState(nf, [1, 3])

        # identity string: <f|c*I|f> = c
        msum = MajoranaSum(Float64, nf)
        PauliPropagation.PropagationBase.add!(msum, getinttype(nf)(0), 0.75)
        @test overlapwithfock(msum, f) == 0.75

        # empty sum
        @test overlapwithfock(MajoranaSum(Float64, nf), f) == 0.0

        # strings with unpaired Majorana modes have zero diagonal overlap
        # (parity superselection), e.g. gamma_1 gamma_3 (both sites unpaired)
        msum = MajoranaSum(Float64, nf)
        PauliPropagation.PropagationBase.add!(msum, MajoranaString(nf, [1, 3]).gammas, 1.0)
        @test overlapwithfock(msum, f) == 0.0

        # fock_mask removes exactly the terms with zero Fock overlap
        msum = random_even_msum(nf, 8)
        masked = fock_mask(msum)
        for fs in allfockstates(nf)
            @test abs(overlapwithfock(msum, fs) - overlapwithfock(masked, fs)) < 1e-14
        end
    end

    @testset "off-diagonal matrix elements <f1|ms|f2>" begin
        nf = 3
        gam = dense_majoranas(nf)
        focks = allfockstates(nf)
        TT = getinttype(nf)
        # arbitrary strings, including odd weight and identity
        strings = [random_string(nf) for _ in 1:25]
        push!(strings, TT(0))
        for ms in strings
            msum = MajoranaSum(ComplexF64, nf)
            PauliPropagation.PropagationBase.add!(msum, ms, 1.0 + 0.0im)
            D = dense_sum(msum, gam)
            for f1 in focks, f2 in focks
                pkg = overlapwithfock(msum, f1, f2)
                dense = dense_fock(f1)' * D * dense_fock(f2)
                @test abs(pkg - dense) < 1e-12
            end
        end
    end

    @testset "spinless superposition overlaps" begin
        nf = 4
        gam = dense_majoranas(nf)
        focks = allfockstates(nf)
        for n_states in 2:4, complexcoeffs in (false, true)
            for _ in 1:10
                msum = random_even_msum(nf, 6)
                D = dense_sum(msum, gam)
                sup_focks = focks[randperm(length(focks))[1:n_states]]
                alphas = complexcoeffs ? randn(ComplexF64, n_states) : complex.(randn(n_states))
                alphas ./= sqrt(sum(abs2, alphas))
                psi = sum(alphas[i] * dense_fock(sup_focks[i]) for i in 1:n_states)
                pkg = overlapwithfock(msum, sup_focks, collect(alphas))
                @test abs(pkg - real(psi' * D * psi)) < 1e-12
            end
        end

        # single-state superposition reduces to the diagonal overlap
        msum = random_even_msum(nf, 6)
        f = focks[7]
        @test abs(overlapwithfock(msum, [f], [1.0 + 0.0im]) - overlapwithfock(msum, f)) < 1e-14

        # unnormalized superposition coefficients are rejected
        @test_throws AssertionError overlapwithfock(msum, focks[1:2], [1.0 + 0.0im, 1.0 + 0.0im])
    end

    @testset "spinful diagonal and superposition" begin
        n_sites = 3
        nf = 2 * n_sites
        gam = dense_majoranas(nf)
        focks = allfockstates(n_sites; is_spinful=true)
        for _ in 1:10
            msum = random_even_msum(nf, 6; is_spinful=true)
            D = dense_sum(msum, gam)
            # diagonal on a sample of Fock states
            for f in focks[randperm(length(focks))[1:12]]
                v = dense_fock(f)
                @test abs(overlapwithfock(msum, f) - real(v' * D * v)) < 1e-12
            end
            # superposition of 3 random states
            sup_focks = focks[randperm(length(focks))[1:3]]
            alphas = randn(ComplexF64, 3)
            alphas ./= sqrt(sum(abs2, alphas))
            psi = sum(alphas[i] * dense_fock(sup_focks[i]) for i in 1:3)
            pkg = overlapwithfock(msum, sup_focks, collect(alphas))
            @test abs(pkg - real(psi' * D * psi)) < 1e-12
        end

        # spinful observables on physical states: <checkerboard|n_up(s)|checkerboard>
        f = FockState(n_sites, [1, 3], [2])   # up on sites 1,3; down on site 2
        for s in 1:n_sites
            @test abs(overlapwithfock(MajoranaSum(n_sites, :nup, s), f) - (s in (1, 3) ? 1.0 : 0.0)) < 1e-14
            @test abs(overlapwithfock(MajoranaSum(n_sites, :ndn, s), f) - (s == 2 ? 1.0 : 0.0)) < 1e-14
        end
    end
end
