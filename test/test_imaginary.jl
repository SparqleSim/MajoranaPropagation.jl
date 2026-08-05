# Imaginary-time evolution (ImaginaryMajoranaRotation / ImaginaryFermionicRotation).
# Conventions under test (docstrings in src/imaginary_gates.jl):
# - A single ImaginaryMajoranaRotation with parameter beta implements
#   rho -> K rho K with K = exp(-beta/2 * G), so on a commuting string
#   ms -> cosh(beta) ms - sinh(beta) G*ms, and anticommuting strings are untouched.
# - ImaginaryFermionicRotation decomposes the symbol operator (identity popped)
#   and applies beta * coeff per Majorana rotation, i.e. K = exp(-beta/2 * (op - c0 I)).
#   Note: unlike real-time FermionicRotation there is NO factor 2 on the angle.
# - normalize_coeffs=true (default) rescales by the identity coefficient after each
#   sub-rotation; this only exists on the ImaginaryFermionicRotation path — a bare
#   ImaginaryMajoranaRotation never normalizes.
using MajoranaPropagation
using PauliPropagation
using LinearAlgebra
using Test
using Random
Random.seed!(42)

if !isdefined(Main, :dense_majoranas)
    include("testutils_dense.jl")
end

@testset "Imaginary-time evolution" begin

    @testset "analytic single-gate action" begin
        nf = 3
        beta = 0.7
        TT = getinttype(nf)
        g = MajoranaString(nf, [1, 2])
        gate = ImaginaryMajoranaRotation(g)

        # evolve the identity: K*I*K = exp(-beta*G) = cosh(beta) I - sinh(beta) G
        rho = MajoranaSum(Float64, nf)
        PauliPropagation.PropagationBase.add!(rho, TT(0), 1.0)
        evolved = propagate!([gate], deepcopy(rho), [beta]; heisenberg=false, min_abs_coeff=-1.0)
        @test length(evolved) == 2
        @test PauliPropagation.PropagationBase.getmergedcoeff(evolved, TT(0)) ≈ cosh(beta)
        @test PauliPropagation.PropagationBase.getmergedcoeff(evolved, g.gammas) ≈ -sinh(beta)

        # thermal expectation of the gate string itself: <G>_beta = -tanh(beta)
        c_id = PauliPropagation.PropagationBase.getmergedcoeff(evolved, TT(0))
        c_g = PauliPropagation.PropagationBase.getmergedcoeff(evolved, g.gammas)
        @test c_g / c_id ≈ -tanh(beta)

        # anticommuting strings are left untouched: gamma_1 gamma_3 anticommutes with G
        ms_anti = MajoranaString(nf, [1, 3])
        @test !MajoranaPropagation.commutes(g, ms_anti)
        rho = MajoranaSum(Float64, nf)
        PauliPropagation.PropagationBase.add!(rho, ms_anti.gammas, 0.8)
        evolved = propagate!([gate], deepcopy(rho), [beta]; heisenberg=false, min_abs_coeff=-1.0)
        @test evolved == rho

        # beta -> 0 is the identity map (up to the sinh(0) = 0 split-off terms,
        # which the default min_abs_coeff truncation removes)
        rho = random_even_msum(nf, 5)
        evolved = propagate!([gate], deepcopy(rho), [0.0]; heisenberg=false)
        @test evolved == rho
    end

    @testset "dense oracle: K rho K" begin
        nf = 3
        gam = dense_majoranas(nf)
        for _ in 1:10
            g = random_even_string(nf)
            rho = random_even_msum(nf, 4)
            beta = 0.3 + rand()
            evolved = propagate!([ImaginaryMajoranaRotation(MajoranaString(nf, Int(g)))], deepcopy(rho), [beta];
                heisenberg=false, min_abs_coeff=-1.0)
            K = exp(-beta / 2 * dense_string(g, nf, gam))
            @test isapprox(dense_sum(evolved, gam), K * dense_sum(rho, gam) * K, atol=1e-11)
        end

        # fermionic wrapper, unnormalized: K = exp(-beta/2 * (op - c0 I))
        for (symb, sites) in [(:nn, [1, 2]), (:n, [2]), (:hop, [1, 2]), (:pair, [2, 3])]
            op = MajoranaSum(nf, symb, sites)
            MajoranaPropagation.pop_id!(op)
            rho = random_even_msum(nf, 3)
            beta = 0.42
            evolved = propagate!([ImaginaryFermionicRotation(symb, sites)], deepcopy(rho), [beta];
                heisenberg=false, normalize_coeffs=false, min_abs_coeff=-1.0)
            K = exp(-beta / 2 * dense_sum(op, gam))
            @test isapprox(dense_sum(evolved, gam), K * dense_sum(rho, gam) * K, atol=1e-11)
        end
    end

    @testset "Fermi-Dirac physics check" begin
        # Imaginary-time evolution of the maximally mixed state under all :n gates
        # prepares the Gibbs state of H = sum_s (n_s - 1/2), whose occupation is
        # the Fermi function <n_s> = 1/(1 + e^beta). The :n strings all commute,
        # so there is no Trotter error and the result is exact.
        nf = 4
        beta = 0.9
        TT = getinttype(nf)
        rho = MajoranaSum(Float64, nf)
        PauliPropagation.PropagationBase.add!(rho, TT(0), 1.0)   # maximally mixed (up to trace)

        circ = [ImaginaryFermionicRotation(:n, [s]) for s in 1:nf]
        betas = fill(beta, nf)
        evolved = propagate!(circ, deepcopy(rho), betas; heisenberg=false, min_abs_coeff=-1.0)

        # with normalize_coeffs=true (default) the identity coefficient is 1,
        # so Tr[rho O]/Tr[rho] = scalarproduct(rho, O)
        @test PauliPropagation.PropagationBase.getmergedcoeff(evolved, TT(0)) ≈ 1.0
        for s in 1:nf
            n_s = MajoranaSum(nf, :n, s)
            @test MajoranaPropagation.scalarproduct(evolved, n_s) ≈ 1 / (1 + exp(beta))
        end
    end

    @testset "dict vs vector backends" begin
        nf = 4
        beta = 0.35
        for _ in 1:5
            rho = random_even_msum(nf, 6)
            circ = [ImaginaryMajoranaRotation(MajoranaString(nf, Int(random_even_string(nf)))) for _ in 1:4]
            betas = fill(beta, length(circ))

            evolved = propagate!(circ, deepcopy(rho), betas; heisenberg=false, min_abs_coeff=-1.0)
            evolved_vec = propagate!(circ, VectorMajoranaSum(deepcopy(rho)), betas; heisenberg=false, min_abs_coeff=-1.0)

            @test length(evolved) == length(evolved_vec)
            @test evolved == evolved_vec
        end
    end

    @testset "error paths" begin
        nf = 3
        g = MajoranaString(nf, [1, 2])
        rho = random_even_msum(nf, 3)

        # imaginary gates are not defined in the Heisenberg picture (the default)
        @test_throws ErrorException propagate!([ImaginaryMajoranaRotation(g)], deepcopy(rho), [0.5])
        @test_throws ErrorException propagate!([ImaginaryFermionicRotation(:n, [1])], deepcopy(rho), [0.5])

        # only even-weight strings can be rotation generators
        @test_throws AssertionError ImaginaryMajoranaRotation(MajoranaString(nf, [1]))
        @test_throws AssertionError ImaginaryMajoranaRotation(MajoranaString(nf, [1, 2, 3]))
    end
end
