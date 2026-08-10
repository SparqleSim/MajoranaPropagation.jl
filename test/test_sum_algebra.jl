# MajoranaSum arithmetic (+, *, scalar *, commutator, scalarproduct, norm)
# against the dense oracle, canonical anticommutation relations expressed
# through the constructors, non-mutation of inputs, and integer-backend sweeps.
using MajoranaPropagation
using PauliPropagation
using LinearAlgebra
using Test
using Random
Random.seed!(42)

if !isdefined(Main, :dense_majoranas)
    include("testutils_dense.jl")
end

@testset "MajoranaSum algebra" begin

    @testset "arithmetic vs dense oracle" begin
        nf = 3
        gam = dense_majoranas(nf)
        for _ in 1:15
            A = random_even_msum(nf, 4)
            B = random_even_msum(nf, 4)
            dA, dB = dense_sum(A, gam), dense_sum(B, gam)

            @test isapprox(dense_sum(A + B, gam), dA + dB, atol=1e-12)
            @test isapprox(dense_sum(A * B, gam), dA * dB, atol=1e-12)
            @test isapprox(dense_sum(2.5 * A, gam), 2.5 * dA, atol=1e-12)
            @test isapprox(dense_sum(MajoranaPropagation.commutator(A, B), gam), (dA * dB - dB * dA), atol=1e-12)

            # scalarproduct(A, B) = Tr[A B]/2^nf, norm(A) = Frobenius norm/2^(nf/2)
            @test isapprox(MajoranaPropagation.scalarproduct(A, B), tr(dA * dB) / 2^nf, atol=1e-12)
            @test isapprox(MajoranaPropagation.norm(A), sqrt(real(tr(dA' * dA)) / 2^nf), atol=1e-12)
        end

        # commuting sums have zero commutator
        A = MajoranaSum(nf, :n, 1)
        B = MajoranaSum(nf, :n, 2)
        @test MajoranaPropagation.norm(MajoranaPropagation.commutator(A, B)) == 0.0
    end

    @testset "canonical anticommutation via constructors" begin
        # {a_i, a†_j} = delta_ij and {a_i, a_j} = 0, expressed with :f/:fdag sums.
        nf = 3
        id_sum = MajoranaSum(Float64, nf)
        PauliPropagation.PropagationBase.add!(id_sum, getinttype(nf)(0), 1.0)
        for i in 1:nf, j in 1:nf
            a_i = MajoranaSum(nf, :f, i)
            adag_j = MajoranaSum(nf, :fdag, j)
            a_j = MajoranaSum(nf, :f, j)

            acomm = a_i * adag_j + adag_j * a_i
            if i == j
                @test MajoranaPropagation.norm(acomm + (-1.0) * id_sum) < 1e-14
            else
                @test MajoranaPropagation.norm(acomm) < 1e-14
            end
            @test MajoranaPropagation.norm(a_i * a_j + a_j * a_i) < 1e-14
        end

        # every Majorana string squares to the identity: S*S = I
        for nf in (3, 10, 40)
            for _ in 1:20
                S = MajoranaSum(Float64, nf)
                PauliPropagation.PropagationBase.add!(S, random_even_string(nf), 1.0)
                sq = S * S
                @test length(sq) == 1
                @test PauliPropagation.PropagationBase.getmergedcoeff(sq, getinttype(nf)(0)) == 1.0
            end
        end
    end

    @testset "coefficient types and demotion" begin
        nf = 3
        A = random_even_msum(nf, 4)

        # products are computed as ComplexF64 (anticommuting strings contribute
        # factors of i even for real inputs) but Hermitian results, e.g. A*A of a
        # real-coefficient sum, get demoted back to Float64
        @test MajoranaPropagation.coefftype(A * A) == Float64

        # a genuinely complex result stays complex: a†a† - aa times a  is not Hermitian
        f = MajoranaSum(nf, :f, 1)
        fdag = MajoranaSum(nf, :fdag, 1)
        @test MajoranaPropagation.coefftype(f) == ComplexF64
        prod = f * fdag
        @test MajoranaPropagation.coefftype(prod) == Float64        # a a† = 1 - n is real
        @test MajoranaPropagation.coefftype(f * MajoranaSum(nf, :n, 2)) == ComplexF64
    end

    @testset "out-of-place ops do not mutate inputs" begin
        nf = 3
        A = random_even_msum(nf, 4)
        B = random_even_msum(nf, 4)
        A0, B0 = deepcopy(A), deepcopy(B)

        A + B
        A * B
        3.0 * A
        MajoranaPropagation.commutator(A, B)
        MajoranaPropagation.scalarproduct(A, B)
        MajoranaPropagation.norm(A)
        @test A == A0
        @test B == B0
    end

    @testset "shape checks and equality" begin
        nf = 3
        A = random_even_msum(nf, 3)
        B = random_even_msum(nf + 1, 3)
        @test_throws ArgumentError A + B
        @test_throws ArgumentError A * B

        # sums on different site counts / spin structure are never equal
        @test MajoranaSum(Float64, 3) != MajoranaSum(Float64, 4)
        @test MajoranaSum(Float64, 4, false) != MajoranaSum(Float64, 4, true)

        # empty and similar
        S = MajoranaPropagation.similar(A)
        @test length(S) == 0
        @test MajoranaPropagation.norm(S) == 0.0
        @test MajoranaPropagation.coefftype(S) == MajoranaPropagation.coefftype(A)
    end

    @testset "real-time propagation vs dense oracle" begin
        nf = 3
        gam = dense_majoranas(nf)

        # Heisenberg action of a MajoranaRotation is U'OU with U = exp(-i theta/2 G).
        # This must FAIL for the opposite convention U O U'.
        for _ in 1:10
            g = random_even_string(nf)
            O = random_even_msum(nf, 4)
            theta = randn()
            evolved = propagate!([MajoranaRotation(MajoranaString(nf, Int(g)))], deepcopy(O), [theta]; min_abs_coeff=-1.0)
            U = exp(-im * theta / 2 * dense_string(g, nf, gam))
            @test isapprox(dense_sum(evolved, gam), U' * dense_sum(O, gam) * U, atol=1e-11)
        end

        # FermionicRotation with angle theta implements exp(-i theta op)
        # (the theta/2 of the underlying MajoranaRotations is cancelled by a factor 2)
        for (symb, sites) in [(:hop, [1, 2]), (:nn, [1, 3]), (:pair, [2, 3])]
            op = MajoranaSum(nf, symb, sites)
            O = random_even_msum(nf, 4)
            theta = 0.37
            evolved = propagate!([FermionicRotation(symb, sites)], deepcopy(O), [theta]; min_abs_coeff=-1.0)
            U = exp(-im * theta * dense_sum(op, gam))
            @test isapprox(dense_sum(evolved, gam), U' * dense_sum(O, gam) * U, atol=1e-11)
        end

        # unitarity: the 2-norm of coefficients is conserved
        O = random_even_msum(nf, 5)
        circ = [MajoranaRotation(MajoranaString(nf, Int(random_even_string(nf)))) for _ in 1:5]
        evolved = propagate!(circ, deepcopy(O), randn(5); min_abs_coeff=-1.0)
        @test isapprox(MajoranaPropagation.norm(evolved), MajoranaPropagation.norm(O), atol=1e-12)

        # round trip: inverse circuit (reversed order, negated angles) restores O
        thetas = randn(5)
        forward = propagate!(circ, deepcopy(O), thetas; min_abs_coeff=-1.0)
        back = propagate!(reverse(circ), forward, -reverse(thetas); min_abs_coeff=-1.0)
        for (ms, coeff) in zip(terms(O), coefficients(O))
            @test isapprox(PauliPropagation.PropagationBase.getmergedcoeff(back, ms), coeff, atol=1e-12)
        end
    end

    @testset "integer backend sweep" begin
        # cross the UInt8/UInt16/.../BitIntegers boundaries and check that the
        # string type is preserved through algebra and propagation
        for nf in (4, 16, 32, 33, 64, 65, 100)
            TT = getinttype(nf)
            msum = random_even_msum(nf, 5)
            @test MajoranaPropagation.majoranatype(msum) == TT

            g = random_even_string(nf)
            @test typeof(g) == TT
            pref, ms3 = ms_mult(MajoranaString(nf, g), MajoranaString(nf, random_even_string(nf)))
            @test typeof(ms3.gammas) == TT
            @test pref in (1, -1, 1im, -1im)

            evolved = propagate!([MajoranaRotation(MajoranaString(nf, g))], deepcopy(msum), [0.3]; min_abs_coeff=-1.0)
            @test MajoranaPropagation.majoranatype(evolved) == TT
            @test MajoranaPropagation.coefftype(evolved) == Float64

            # dict and vector backends agree at every width
            evolved_vec = propagate!([MajoranaRotation(MajoranaString(nf, g))], VectorMajoranaSum(deepcopy(msum)), [0.3]; min_abs_coeff=-1.0)
            @test evolved == evolved_vec
        end
    end
end
