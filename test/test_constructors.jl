# Symbol constructors (src/Constructors.jl) pinned against independently built
# second-quantized dense operators. These tests define the physical meaning of
# every symbol; if one fails, either the constructor or the documented
# convention has drifted.
using MajoranaPropagation
using PauliPropagation
using LinearAlgebra
using Test
using Random
Random.seed!(42)

if !isdefined(Main, :dense_majoranas)
    include("testutils_dense.jl")
end

@testset "Symbol constructors vs second quantization" begin

    @testset "spinless operators" begin
        nf = 3
        gam = dense_majoranas(nf)
        a = dense_annihilators(nf, gam)
        ad = [Matrix(x') for x in a]
        n = [ad[s] * a[s] for s in 1:nf]

        for s in 1:nf
            @test isapprox(dense_sum(MajoranaSum(nf, :n, s), gam), n[s], atol=1e-13)
            @test isapprox(dense_sum(MajoranaSum(nf, :f, s), gam), a[s], atol=1e-13)
            @test isapprox(dense_sum(MajoranaSum(nf, :fdag, s), gam), ad[s], atol=1e-13)
        end
        for (i, j) in [(1, 2), (1, 3), (2, 3)]
            @test isapprox(dense_sum(MajoranaSum(nf, :hop, [i, j]), gam), ad[i] * a[j] + ad[j] * a[i], atol=1e-13)
            @test isapprox(dense_sum(MajoranaSum(nf, :nn, [i, j]), gam), n[i] * n[j], atol=1e-13)
            # pair creation convention: a†_i a†_j - a_i a_j (i < j)
            @test isapprox(dense_sum(MajoranaSum(nf, :pair, [i, j]), gam), ad[i] * ad[j] - a[i] * a[j], atol=1e-13)
        end
    end

    @testset "spinful operators" begin
        n_sites = 2
        nf = 2 * n_sites
        gam = dense_majoranas(nf)
        a = dense_annihilators(nf, gam)
        ad = [Matrix(x') for x in a]
        # mode layout: site s -> up at 2s-1, down at 2s
        up(s) = 2s - 1
        dn(s) = 2s
        n_up = [ad[up(s)] * a[up(s)] for s in 1:n_sites]
        n_dn = [ad[dn(s)] * a[dn(s)] for s in 1:n_sites]
        id = Matrix{ComplexF64}(I, 2^nf, 2^nf)

        for s in 1:n_sites
            @test isapprox(dense_sum(MajoranaSum(n_sites, :nup, s), gam), n_up[s], atol=1e-13)
            @test isapprox(dense_sum(MajoranaSum(n_sites, :ndn, s), gam), n_dn[s], atol=1e-13)
            @test isapprox(dense_sum(MajoranaSum(n_sites, :nupndn, s), gam), n_up[s] * n_dn[s], atol=1e-13)
            @test isapprox(dense_sum(MajoranaSum(n_sites, :hole, s), gam), (id - n_up[s]) * (id - n_dn[s]), atol=1e-13)
            @test isapprox(dense_sum(MajoranaSum(n_sites, :Sz, s), gam), (n_up[s] - n_dn[s]) / 2, atol=1e-13)
            # on-site spin flip
            @test isapprox(dense_sum(MajoranaSum(n_sites, :hop_on_site, s), gam),
                ad[up(s)] * a[dn(s)] + ad[dn(s)] * a[up(s)], atol=1e-13)
            @test isapprox(dense_sum(MajoranaSum(n_sites, :fup, s), gam), a[up(s)], atol=1e-13)
            @test isapprox(dense_sum(MajoranaSum(n_sites, :fupdag, s), gam), ad[up(s)], atol=1e-13)
            @test isapprox(dense_sum(MajoranaSum(n_sites, :fdn, s), gam), a[dn(s)], atol=1e-13)
            @test isapprox(dense_sum(MajoranaSum(n_sites, :fdndag, s), gam), ad[dn(s)], atol=1e-13)
        end

        i, j = 1, 2
        @test isapprox(dense_sum(MajoranaSum(n_sites, :hopup, [i, j]), gam),
            ad[up(i)] * a[up(j)] + ad[up(j)] * a[up(i)], atol=1e-13)
        @test isapprox(dense_sum(MajoranaSum(n_sites, :hopdn, [i, j]), gam),
            ad[dn(i)] * a[dn(j)] + ad[dn(j)] * a[dn(i)], atol=1e-13)
        @test isapprox(dense_sum(MajoranaSum(n_sites, :pairup, [i, j]), gam),
            ad[up(i)] * ad[up(j)] - a[up(i)] * a[up(j)], atol=1e-13)
        @test isapprox(dense_sum(MajoranaSum(n_sites, :pairdn, [i, j]), gam),
            ad[dn(i)] * ad[dn(j)] - a[dn(i)] * a[dn(j)], atol=1e-13)
        # :hopupdn couples (site i, up) to (site j, down) and must NOT reorder sites
        @test isapprox(dense_sum(MajoranaSum(n_sites, :hopupdn, [i, j]), gam),
            ad[up(i)] * a[dn(j)] + ad[dn(j)] * a[up(i)], atol=1e-13)
        # BUG: when the up-site index is larger than the down-site index the
        # constructor currently returns MINUS the hopping operator (the stored
        # +-0.5 pattern assumes ascending Majorana mode order, which flips here).
        # Flip this to @test once fixed in src/Constructors.jl.
        @test_broken isapprox(dense_sum(MajoranaSum(n_sites, :hopupdn, [j, i]), gam),
            ad[up(j)] * a[dn(i)] + ad[dn(i)] * a[up(j)], atol=1e-13)
    end

    @testset "site ordering" begin
        nf = 3
        # two-site symbols sort their site indices ...
        for symb in (:hop, :nn, :pair)
            @test MajoranaSum(nf, symb, [3, 1]) == MajoranaSum(nf, symb, [1, 3])
        end
        # ... except :hopupdn, whose indices refer to different spins
        @test MajoranaSum(2, :hopupdn, [2, 1]) != MajoranaSum(2, :hopupdn, [1, 2])
    end

    @testset "physics invariants" begin
        nf = 3
        # total particle number N = sum_s n_s
        N = MajoranaSum(nf, :n, 1) + MajoranaSum(nf, :n, 2) + MajoranaSum(nf, :n, 3)
        # hopping and density terms conserve particle number, pairing does not
        @test MajoranaPropagation.norm(MajoranaPropagation.commutator(MajoranaSum(nf, :hop, [1, 2]), N)) < 1e-14
        @test MajoranaPropagation.norm(MajoranaPropagation.commutator(MajoranaSum(nf, :nn, [1, 2]), N)) < 1e-14
        @test MajoranaPropagation.norm(MajoranaPropagation.commutator(MajoranaSum(nf, :pair, [1, 2]), N)) > 0.1

        # all observable-type constructors have real coefficients (Hermitian operators)
        for (symb, sites) in [(:n, 1), (:hop, [1, 2]), (:nn, [1, 2]), (:pair, [1, 2])]
            @test MajoranaPropagation.coefftype(MajoranaSum(nf, symb, sites)) == Float64
        end
    end

    @testset "fermionic gate decomposition" begin
        n_sites = 4
        # :nn = 4 stored terms, identity popped -> 3 rotations, truncate after each
        rotations, coeffs, truncate_each = MajoranaPropagation.getmajoranarotations(FermionicRotation(:nn, [1, 2]), n_sites)
        @test length(rotations) == 3
        @test length(coeffs) == 3
        @test truncate_each

        # :hop = 2 terms, no identity -> 2 rotations, flagged non-number-preserving
        rotations, coeffs, truncate_each = MajoranaPropagation.getmajoranarotations(FermionicRotation(:hop, [1, 2]), n_sites)
        @test length(rotations) == 2
        @test sort(coeffs) == [-0.5, 0.5]
        @test !truncate_each

        @test MajoranaPropagation.flag_non_number_preserving(:hop)
        @test MajoranaPropagation.flag_non_number_preserving(:hopup)
        @test MajoranaPropagation.flag_non_number_preserving(:hopdn)
        @test !MajoranaPropagation.flag_non_number_preserving(:nn)
        @test !MajoranaPropagation.flag_non_number_preserving(:n)

        # convenience constructors
        @test FermionicRotation(:n, 2).sites == [2]
        @test FermionicRotation(:hop, (1, 3)).sites == [1, 3]
    end

    @testset "error paths" begin
        @test_throws ErrorException MajoranaSum(3, :bogus, 1)
        @test_throws AssertionError MajoranaSum(3, :hop, [1, 2, 3])
    end
end
