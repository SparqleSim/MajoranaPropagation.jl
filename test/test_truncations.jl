# Truncation logic: unit tests of the predicates, soundness (loose thresholds
# must not change results), monotonicity, and a pinned known-value regression.
using MajoranaPropagation
using PauliPropagation
using Test
using Random
Random.seed!(42)

if !isdefined(Main, :dense_majoranas)
    include("testutils_dense.jl")
end

# deterministic spinful Hubbard Trotter circuit used in several testsets
function _hubbard_circ(N_sites; t=1.0, U=1.0, dtau=0.1)
    topo = bricklayertopology(N_sites)
    circ = FermionicRotation[]
    thetas = Float64[]
    for pair in topo
        push!(circ, FermionicRotation(:hopup, pair))
        push!(thetas, -t * dtau)
    end
    for pair in topo
        push!(circ, FermionicRotation(:hopdn, pair))
        push!(thetas, -t * dtau)
    end
    for i in 1:N_sites
        push!(circ, FermionicRotation(:nupndn, i))
        push!(thetas, U * dtau)
    end
    return circ, thetas
end

@testset "Truncations" begin

    @testset "weight and unpaired predicates" begin
        nf = 4
        mask = create_unpaired_mask(nf)

        paired = MajoranaString(nf, [1, 2]).gammas            # gamma_1 gamma'_1
        unpaired2 = MajoranaString(nf, [1, 3]).gammas         # gamma_1 gamma_2
        mixed = MajoranaString(nf, [1, 2, 3]).gammas          # site 1 paired, site 2 unpaired
        idstr = getinttype(nf)(0)

        @test compute_unpaired(paired, mask) == 0
        @test compute_unpaired(unpaired2, mask) == 2
        @test compute_unpaired(mixed, mask) == 1
        @test compute_unpaired(idstr, mask) == 0

        @test MajoranaPropagation.truncateunpaired(unpaired2, 1, mask)
        @test !MajoranaPropagation.truncateunpaired(unpaired2, 2, mask)

        w4 = MajoranaString(nf, [1, 2, 3, 4]).gammas
        @test truncatemajoranaweight(w4, 3)
        @test !truncatemajoranaweight(w4, 4)
        @test !truncatemajoranaweight(idstr, 0)
    end

    @testset "unpaired mask across integer widths" begin
        # getinttype spans UInt8 .. UInt128 and BitIntegers types beyond;
        # the mask must be exactly the odd bit positions at every width
        for nf in (1, 2, 3, 4, 8, 16, 31, 32, 33, 64, 65)
            TT = getinttype(nf)
            mask = create_unpaired_mask(nf)
            @test mask isa TT
            @test mask == create_unpaired_mask(TT, nf)
            @test mask == reduce(|, (TT(1) << (2k + 1) for k in 0:nf-1))

            # compute_unpaired agrees with a per-mode bit count on random strings
            for _ in 1:5
                ms = TT(0)
                for i in 0:2nf-1
                    ms |= TT(rand(Bool)) << i
                end
                expected = count(k -> count_ones((ms >>> (2k)) & TT(3)) == 1, 0:nf-1)
                @test compute_unpaired(ms, mask) == expected
            end
        end
    end

    @testset "doublon counters (spinful)" begin
        n_sites = 3
        filters = create_doublons_filters(n_sites)
        nf = 2 * n_sites

        # all four Majoranas on site 2: indices 4*2-3 .. 4*2
        doublon2 = MajoranaString(nf, [5, 6, 7, 8]).gammas
        @test compute_doublons(doublon2, filters) == 1

        # sites 1 and 3 fully occupied
        doublon13 = MajoranaString(nf, [1, 2, 3, 4, 9, 10, 11, 12]).gammas
        @test compute_doublons(doublon13, filters) == 2

        # partial site occupation is not a doublon
        partial = MajoranaString(nf, [5, 6, 7]).gammas
        @test compute_doublons(partial, filters) == 0
        @test compute_doublons(getinttype(nf)(0), filters) == 0
    end

    @testset "truncate! on a MajoranaSum" begin
        nf = 4
        TT = getinttype(nf)
        msum = MajoranaSum(Float64, nf)
        PauliPropagation.PropagationBase.add!(msum, MajoranaString(nf, [1, 2]).gammas, 1.0)          # w2 paired
        PauliPropagation.PropagationBase.add!(msum, MajoranaString(nf, [1, 3]).gammas, 0.5)          # w2 unpaired
        PauliPropagation.PropagationBase.add!(msum, MajoranaString(nf, [1, 2, 3, 4]).gammas, 1e-6)   # w4 small coeff
        PauliPropagation.PropagationBase.add!(msum, MajoranaString(nf, [1, 2, 3, 4, 5, 6]).gammas, 2.0) # w6

        # min_abs_coeff removes exactly the small term
        trunc = truncate!(deepcopy(msum); min_abs_coeff=1e-3)
        @test length(trunc) == 3
        @test PauliPropagation.PropagationBase.getmergedcoeff(trunc, MajoranaString(nf, [1, 2, 3, 4]).gammas) == 0.0

        # max_weight removes the weight-6 (and, with min_abs_coeff, the small) term
        trunc = truncate!(deepcopy(msum); min_abs_coeff=-1.0, max_weight=4)
        @test length(trunc) == 3
        @test PauliPropagation.PropagationBase.getmergedcoeff(trunc, MajoranaString(nf, [1, 2, 3, 4, 5, 6]).gammas) == 0.0

        # max_unpaired=0 keeps exactly the terms fock_mask keeps
        trunc = truncate!(deepcopy(msum); min_abs_coeff=-1.0, max_unpaired=0)
        masked = fock_mask(msum)
        @test sort(collect(terms(trunc))) == sort(collect(terms(masked)))

        # customtruncfunc: drop everything overlapping site 1 (bits 1 and 2)
        msum2 = MajoranaSum(Float64, nf)
        PauliPropagation.PropagationBase.add!(msum2, MajoranaString(nf, [1, 2]).gammas, 1.0)   # touches site 1
        PauliPropagation.PropagationBase.add!(msum2, MajoranaString(nf, [2, 5]).gammas, 1.0)   # touches site 1
        PauliPropagation.PropagationBase.add!(msum2, MajoranaString(nf, [3, 4]).gammas, 1.0)   # site 2 only
        site1mask = MajoranaString(nf, [1, 2]).gammas
        trunc = truncate!(deepcopy(msum2); min_abs_coeff=-1.0, customtruncfunc=(mstr, coeff) -> (mstr & site1mask) != 0)
        @test length(trunc) == 1
        @test PauliPropagation.PropagationBase.getmergedcoeff(trunc, MajoranaString(nf, [3, 4]).gammas) == 1.0
    end

    @testset "loose thresholds do not change results" begin
        N_sites = 4
        nf = 2 * N_sites
        circ, thetas = _hubbard_circ(N_sites)
        obs0 = MajoranaSum(N_sites, :nupndn, 2)
        MajoranaPropagation.pop_id!(obs0)

        baseline = propagate!(circ, deepcopy(obs0), thetas; min_abs_coeff=-1.0)
        loose = propagate!(circ, deepcopy(obs0), thetas;
            min_abs_coeff=-1.0, max_weight=2nf, max_unpaired=2nf, max_freq=Inf, max_sins=Inf)
        @test baseline == loose
    end

    @testset "monotonicity in truncation strength" begin
        N_sites = 4
        circ, thetas = _hubbard_circ(N_sites)
        fock = FockState(N_sites, :checkerboard, true)
        obs0 = MajoranaSum(N_sites, :nupndn, 2)
        MajoranaPropagation.pop_id!(obs0)

        baseline = propagate!(circ, deepcopy(obs0), thetas; min_abs_coeff=-1.0)
        exact_val = overlapwithfock(baseline, fock)

        # tightening min_abs_coeff never increases the term count,
        # and the error vs the untruncated result shrinks as the threshold loosens
        cutoffs = [1e-1, 1e-2, 1e-3, 1e-4]
        lengths = Int[]
        errors = Float64[]
        for cutoff in cutoffs
            obs = propagate!(circ, deepcopy(obs0), thetas; min_abs_coeff=cutoff)
            push!(lengths, length(obs))
            push!(errors, abs(overlapwithfock(obs, fock) - exact_val))
        end
        @test issorted(lengths)
        @test issorted(errors, rev=true) || all(errors .< 1e-12)

        # same for max_weight
        lengths = Int[]
        for w in (2, 4, 6, 8)
            obs = propagate!(circ, deepcopy(obs0), thetas; min_abs_coeff=-1.0, max_weight=w)
            push!(lengths, length(obs))
        end
        @test issorted(lengths)
        # and evolved sums respect the weight cap
        obs = propagate!(circ, deepcopy(obs0), thetas; min_abs_coeff=-1.0, max_weight=4)
        @test all(get_weight(ms) <= 4 for ms in terms(obs))
    end

    @testset "pinned known-value regression" begin
        # Deterministic circuit and parameters; the value below was cross-checked
        # against the dense/Yao-verified paths when this test was written.
        # If this fails after a code change, a physics convention has drifted.
        N_sites = 6
        circ, thetas = _hubbard_circ(N_sites)
        obs = MajoranaSum(N_sites, :nupndn, 3)
        MajoranaPropagation.pop_id!(obs)
        fock = FockState(N_sites, :checkerboard, true)
        for _ in 1:4
            propagate!(circ, obs, thetas; min_abs_coeff=1e-5)
        end
        expected_value = -0.054894655344715965
        @test isapprox(overlapwithfock(obs, fock), expected_value; rtol=1e-6)
    end

    # Out-of-place `propagate` vs in-place `propagate!` (TESTING_GUIDE.md §4):
    # equal results with ==, no mutation of the input, on both backends.
    # Regression test for the UndefVarError in `propagate` introduced in the
    # performance-updates PR (#22), whose @test_broken marker was removed there.
    @testset "propagate vs propagate! ($label)" for (label, make) in
        (("dict", identity), ("vector", VectorMajoranaSum))

        N_sites = 3
        circ, thetas = _hubbard_circ(N_sites)
        msum0 = make(MajoranaSum(N_sites, :nupndn, 2))
        msum_in = deepcopy(msum0)
        res = propagate(circ, msum_in, thetas; min_abs_coeff=1e-12)
        @test msum_in == msum0
        @test res == propagate!(circ, deepcopy(msum0), thetas; min_abs_coeff=1e-12)

        # max_freq/max_sins on plain-number coefficients: propagate must auto-wrap into
        # MajoranaFrequencyTracker, truncate, and hand back unwrapped plain numbers
        res = propagate(circ, msum_in, thetas; min_abs_coeff=-1.0, max_freq=3)
        @test msum_in == msum0
        @test MajoranaPropagation.coefftype(res) == Float64

        ref = propagate!(circ, wrapcoefficients(deepcopy(msum0), MajoranaFrequencyTracker),
            thetas; min_abs_coeff=-1.0, max_freq=3)
        @test res == MajoranaPropagation.unwrapcoefficients(ref)

        # the truncation must actually bite compared to the untruncated run
        @test length(res) < length(propagate(circ, msum0, thetas; min_abs_coeff=-1.0))
    end
end
