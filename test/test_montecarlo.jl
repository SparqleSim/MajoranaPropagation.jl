# Monte Carlo propagation: `mcsample` keeps one random branch per Majorana string and rotation, `mcpropagate`
# resamples the sum down past `max_size`. The loops and the resampling strategies are PropagationBase's; here
# the sampling rule of a Majorana rotation is checked against the deterministic branches, and every backend
# runs both.
using MajoranaPropagation
using PauliPropagation
using Test
using Random

if !isdefined(Main, :dense_majoranas)
    include("testutils_dense.jl")
end
# `hubbard_layer`, `backends`, `termdict` and `agree` come from test_backends.jl
if !isdefined(Main, :hubbard_layer)
    include("test_backends.jl")
end

@testset "Monte Carlo" begin

    @testset "mcapplytoall! keeps one of the two exact branches" begin
        nf = 2
        gate_ms = MajoranaString(nf, [1, 2])
        gate = MajoranaRotation(gate_ms)
        term = MajoranaString(nf, [1, 3]).gammas
        @test !MajoranaPropagation.commutes(gate_ms.gammas, term)
        gps = MajoranaPropagation.compute_parity_bits_and_shift(gate_ms.gammas, 2nf)
        new_term, prod_sign = majoranarotationproduct(term, gate_ms.gammas, gps)
        theta = 0.37
        cos_val, sin_val = cos(theta), sin(theta)
        normalization = abs(cos_val) + abs(sin_val)

        # one call keeps one branch, boosted by 1/probability so that averaging many calls gives the split
        seen_stay = seen_flip = false
        for _ in 1:200
            vmsum = VectorMajoranaSum(nf, false, [term], [1.0])
            mcapplytoall!(gate, vmsum, theta)
            t, c = only(terms(vmsum)), only(coefficients(vmsum))
            if t == term
                @test c ≈ normalization * sign(cos_val)
                seen_stay = true
            else
                @test t == new_term
                @test c ≈ normalization * sign(sin_val) * prod_sign
                seen_flip = true
            end
        end
        @test seen_stay && seen_flip

        # a commuting string is left untouched
        commuting = MajoranaString(nf, [3, 4]).gammas
        vmsum = VectorMajoranaSum(nf, false, [commuting], [1.0])
        mcapplytoall!(gate, vmsum, theta)
        @test only(terms(vmsum)) == commuting && only(coefficients(vmsum)) == 1.0

        # squared=true samples by cos^2 and sin^2, which add to 1, so |coeff| is preserved exactly
        seen_stay = seen_flip = false
        for _ in 1:200
            vmsum = VectorMajoranaSum(nf, false, [term], [2.5])
            mcapplytoall!(gate, vmsum, 0.9; squared=true)
            @test only(coefficients(vmsum)) == 2.5
            seen_stay |= only(terms(vmsum)) == term
            seen_flip |= only(terms(vmsum)) == new_term
        end
        @test seen_stay && seen_flip
    end

    @testset "every backend keeps one branch per string" begin
        nf = 3
        gate_ms = MajoranaString(nf, [1, 2, 3, 4])
        gate = MajoranaRotation(gate_ms)
        gps = MajoranaPropagation.compute_parity_bits_and_shift(gate_ms.gammas, 2nf)
        # three anticommuting strings and one commuting, whose branches are all distinct, so no backend merges
        strings = [MajoranaString(nf, v).gammas for v in ([1, 5], [3, 6], [1, 2, 3, 5], [5, 6])]
        branches = Set(vcat(strings, [first(majoranarotationproduct(s, gate_ms.gammas, gps)) for s in strings]))
        @test length(branches) == 8
        obs = MajoranaSum(nf, false, Dict(s => 1.0 for s in strings))

        for b in backends(obs)
            sampled = mcapplytoall!(gate, deepcopy(b), 0.7; squared=true, thread=false)
            @test length(sampled) == 4
            @test all(t in branches for (t, _) in sampled)
            @test all(abs(c) ≈ 1 for c in coefficients(sampled))
        end

        # path properties count the branch a string keeps
        Random.seed!(7)
        plain = mcapplytoall!(gate, VectorMajoranaSum(obs), 0.7; thread=false)
        Random.seed!(7)
        tracked = mcapplytoall!(gate, wrapcoefficients(VectorMajoranaSum(obs), MajoranaFrequencyTracker), 0.7; thread=false)
        @test terms(tracked) == terms(plain)
        @test [c.coeff for c in coefficients(tracked)] ≈ coefficients(plain)
        @test all(c.freq == c.nsins + c.ncos == !MajoranaPropagation.commutes(gate_ms.gammas, s)
                  for (s, c) in zip(terms(VectorMajoranaSum(obs)), coefficients(tracked)))
    end

    @testset "mcsample statistically reproduces propagate" begin
        Random.seed!(1)
        N = 4
        circ, thetas = hubbard_layer(bricklayertopology(N), N)
        obs = MajoranaSum(N, :nupndn, 2)
        fock = FockState(N, :checkerboard, true)
        exact = overlapwithfock(propagate(circ, obs, thetas; min_abs_coeff=0.0), fock)

        # `reps` walkers per string of the observable, each sampling its own path; the sum over all walkers
        # divided by `reps` is the unbiased estimate
        reps = 20_000
        strings, coeffs = collect(terms(obs)), collect(coefficients(obs))
        init = VectorMajoranaSum(N, true, repeat(strings, reps), repeat(coeffs, reps))
        # a multi sum merges what it is built from, so its zones are replicated one by one
        base = MultiMajoranaSum(VectorMajoranaSum(obs), 4)
        multi = MultiMajoranaSum(N, [VectorMajoranaSum(N, true, repeat(z.terms, reps), repeat(z.coeffs, reps)) for z in zones(base)], MajoranaPropagation.zonemap(base))
        @test length(multi) == length(init)
        for msum in (init, multi)
            sampled = mcsample(circ, msum, thetas)
            @test length(sampled) == length(init)
            @test overlapwithfock(sampled, fock) / reps ≈ exact atol = 0.05
        end
    end

    @testset "mcpropagate below max_size is propagate, on every backend" begin
        N = 6
        circ, thetas = hubbard_layer(bricklayertopology(N), N)
        circ, thetas = repeat(circ, 2), repeat(thetas, 2)
        obs = MajoranaSum(N, :nupndn, 3)
        min_abs_coeff = 1e-5

        # max_size never hit: applymergetruncateresample! is applymergetruncate!, bit for bit
        exact = propagate(circ, obs, thetas; min_abs_coeff)
        results = [mcpropagate(circ, b, thetas; max_size=10^9, min_abs_coeff) for b in backends(obs)]
        @test agree([exact; results])
        @test [typeof(r) for r in results] == [typeof(b) for b in backends(obs)]
        @test length(exact) > 10

        # past max_size the sum is resampled down, through either default strategy
        for b in backends(obs), squared in (false, true)
            result = mcpropagate(circ, b, thetas; max_size=10, squared, min_abs_coeff=1e-8)
            @test 0 < length(result) <= 10
        end
    end

    @testset "mcpropagate through imaginary-time gates" begin
        N = 4
        circ = ImaginaryFermionicRotation[]
        taus = Float64[]
        for pair in bricklayertopology(N)
            push!(circ, ImaginaryFermionicRotation(:hopup, pair))
            push!(taus, -0.05)
            push!(circ, ImaginaryFermionicRotation(:hopdn, pair))
            push!(taus, -0.05)
        end
        for i in 1:N
            push!(circ, ImaginaryFermionicRotation(:nupndn, [i]))
            push!(taus, 0.2)
        end
        circ, taus = repeat(circ, 3), repeat(taus, 3)
        rho0 = MajoranaSum(Float64, N, true)
        add!(rho0, getinttype(2N)(0), 1.0)

        exact = propagate(circ, rho0, taus; heisenberg=false, min_abs_coeff=1e-4)
        results = [mcpropagate(circ, b, taus; max_size=10^9, heisenberg=false, min_abs_coeff=1e-4) for b in backends(rho0)]
        @test agree([exact; results])
        @test length(exact) > 10

        # resampling may drop the identity, so the bounded run does not normalize by it
        for b in backends(rho0)
            result = mcpropagate(circ, b, taus; max_size=10, heisenberg=false, normalize_coeffs=false, min_abs_coeff=1e-4)
            @test 0 < length(result) <= 10
        end

        # path sampling has no rule for imaginary-time gates
        @test_throws ErrorException mcsample(circ, rho0, taus; heisenberg=false)
    end

    @testset "the circuit is applied in the Heisenberg picture" begin
        # two anticommuting rotations by pi/2: each flips a string that anticommutes with it for certain, after
        # which the string commutes with the other, so the gate the Heisenberg order applies first decides the result
        nf = 2
        g1 = MajoranaRotation(MajoranaString(nf, [1, 2]))
        g2 = MajoranaRotation(MajoranaString(nf, [2, 3]))
        obs = MajoranaSum(nf, [1, 3], false)
        flipped_by_g2 = MajoranaString(nf, [1, 2]).gammas
        flipped_by_g1 = MajoranaString(nf, [2, 3]).gammas

        for b in backends(obs)
            sampled = mcsample([g1, g2], b, [pi / 2, pi / 2])
            @test termdict(sampled) == termdict(propagate([g1, g2], b, [pi / 2, pi / 2]))
            @test collect(keys(termdict(sampled))) == [flipped_by_g2]
            @test collect(keys(termdict(mcsample([g2, g1], b, [pi / 2, pi / 2])))) == [flipped_by_g1]
        end

        # a single gate and parameter are promoted to lists, and the Schrödinger picture is not defined
        @test termdict(mcsample(g1, obs, pi / 2)) == termdict(propagate(g1, obs, pi / 2))
        @test_throws ErrorException mcsample(g1, obs, 0.3; heisenberg=false)
        @test_throws ErrorException mcpropagate(g1, obs, 0.3; max_size=10, heisenberg=false)
    end

    @testset "in place and out of place" begin
        nf = 3
        gate = MajoranaRotation(MajoranaString(nf, [1, 2]))
        obs = random_even_msum(nf, 5)
        for b in backends(obs)
            original = deepcopy(b)
            @test mcsample!(gate, b, 0.3) === b
            fresh = deepcopy(original)
            @test mcsample(gate, fresh, 0.3) !== fresh
            @test termdict(fresh) == termdict(original)
            @test mcpropagate!(gate, b, 0.3; max_size=100) === b
            fresh = deepcopy(original)
            @test mcpropagate(gate, fresh, 0.3; max_size=100) !== fresh
            @test termdict(fresh) == termdict(original)
        end

        # on a vector or dict cache, mcsample! samples the active sum in place and never swaps the sums
        for cache in (PropagationCache(VectorMajoranaSum(obs)), PropagationCache(deepcopy(obs)))
            main_before, aux_before = mainsum(cache), auxsum(cache)
            @test mcsample!(gate, cache, 0.3) === cache
            @test mainsum(cache) === main_before && auxsum(cache) === aux_before
            @test 0 < length(cache) <= length(obs)
        end
        @test_throws ArgumentError mcsample!(gate, PropagationCache(MultiMajoranaSum(VectorMajoranaSum(obs), 4)), 0.3)
    end
end
