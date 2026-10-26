# Every backend of a Majorana sum propagates by the same rules through PropagationBase's kernels, so a
# MajoranaSum, a VectorMajoranaSum and a MultiMajoranaSum over either must end with the same Majorana strings
# and the same coefficients, exactly, under every truncation. Run the suite with JULIA_NUM_THREADS>1 to
# exercise the threaded kernels.
using MajoranaPropagation
using PauliPropagation
using LinearAlgebra
using Test
using Random

if !isdefined(Main, :dense_majoranas)
    include("testutils_dense.jl")
end

const backend_rng = MersenneTwister(7)

# one first-order Trotter layer of the spinful Hubbard model: hoppings decompose into two rotations
# that are only number-preserving together, `:nupndn` into rotations truncated one by one
function hubbard_layer(topo, n_sites; t=1.0, U=1.0, dt=0.1)
    circ = FermionicRotation[]
    thetas = Float64[]
    for pair in topo
        push!(circ, FermionicRotation(:hopup, pair))
        push!(thetas, -t * dt)
        push!(circ, FermionicRotation(:hopdn, pair))
        push!(thetas, -t * dt)
    end
    for i in 1:n_sites
        push!(circ, FermionicRotation(:nupndn, i))
        push!(thetas, U * dt)
    end
    return circ, thetas
end

# random even-weight Majorana rotations of type `R` on `nf` spinless fermions
function random_rotation_circuit(nf, n_gates; R=MajoranaRotation)
    circ = R{getinttype(nf)}[]
    thetas = Float64[]
    for _ in 1:n_gates
        w = rand(backend_rng, (2, 4))
        gammas = Int[]
        while length(gammas) < w
            g = rand(backend_rng, 1:2*nf)
            g in gammas || push!(gammas, g)
        end
        push!(circ, R(MajoranaString(nf, gammas)))
        push!(thetas, randn(backend_rng))
    end
    return circ, thetas
end

# O(N) exact comparison of the strings and coefficients of a sum of any backend (the generic `==` re-merges
# copies and looks the strings up one by one)
termdict(msum) = Dict(k => v for (k, v) in msum)

# the backends of the same sum; the multi sums have more zones than a test usually has threads, so that strings
# really move between zones
backends(msum::MajoranaSum; n_zones=4) =
    (msum, VectorMajoranaSum(msum), MultiMajoranaSum(VectorMajoranaSum(msum), n_zones), MultiMajoranaSum(msum, n_zones))

# propagate `obs` from every backend and return the results, the dict one first
propagateall(circ, obs, thetas; kwargs...) = [propagate(circ, b, thetas; kwargs...) for b in backends(obs)]

agree(results) = all(termdict(r) == termdict(first(results)) for r in results)

@testset "Backends" begin

    @testset "every backend agrees under every truncation" begin
        N = 6
        circ, thetas = hubbard_layer(bricklayertopology(N), N)
        circ, thetas = repeat(circ, 2), repeat(thetas, 2)
        obs = MajoranaSum(N, :nupndn, 3)

        for min_abs_coeff in (0.0, 1e-4), max_weight in (Inf, 4, 6), max_unpaired in (Inf, 2, 4)
            results = propagateall(circ, obs, thetas; min_abs_coeff, max_weight, max_unpaired)
            @test agree(results)
            @test length(first(results)) > 1
        end

        # both ways of truncating a fermionic gate's rotations
        for truncate_each_mr in (true, false)
            @test agree(propagateall(circ, obs, thetas; min_abs_coeff=1e-4, max_unpaired=4, truncate_each_mr))
        end

        # a precomputed unpaired mask and a custom truncation
        let unpaired_mask = create_unpaired_mask(2N), customtruncfunc = (mstr, coeff) -> get_weight(mstr) == 4 && abs(coeff) < 1e-2
            @test agree(propagateall(circ, obs, thetas; min_abs_coeff=1e-5, max_unpaired=4, unpaired_mask, customtruncfunc))
        end

        # results come back in the backend they went in
        @test [typeof(r) for r in propagateall(circ, obs, thetas)] == [typeof(b) for b in backends(obs)]
    end

    @testset "every backend agrees for bare Majorana rotations" begin
        nf = 8
        circ, thetas = random_rotation_circuit(nf, 40)
        obs = random_even_msum(nf, 5)
        truncate!(obs; min_abs_coeff=1e-6, max_weight=4)

        for max_weight in (Inf, 4), min_abs_coeff in (0.0, 1e-6)
            @test agree(propagateall(circ, obs, thetas; min_abs_coeff, max_weight))
        end
    end

    @testset "strings dropped as they are made are those the truncation drops" begin
        # `max_weight` and `max_unpaired` keep a rotation from making the strings they reject, while a
        # `customtruncfunc` with the same test drops them only after they are merged
        N = 6
        circ, thetas = hubbard_layer(bricklayertopology(N), N)
        circ, thetas = repeat(circ, 2), repeat(thetas, 2)
        obs = MajoranaSum(N, :nupndn, 3)
        mask = create_unpaired_mask(2N)

        for (kw, customtruncfunc) in (((; max_unpaired=2), (mstr, coeff) -> compute_unpaired(mstr, mask) > 2),
                ((; max_weight=4), (mstr, coeff) -> get_weight(mstr) > 4)), truncate_each_mr in (nothing, true)
            for b in backends(obs)
                dropped = propagate(circ, b, thetas; min_abs_coeff=1e-5, truncate_each_mr, kw...)
                truncated = propagate(circ, b, thetas; min_abs_coeff=1e-5, truncate_each_mr, customtruncfunc)
                @test termdict(dropped) == termdict(truncated)
            end
        end

        nf = 8
        rcirc, rthetas = random_rotation_circuit(nf, 40)
        robs = random_even_msum(nf, 5)
        truncate!(robs; min_abs_coeff=1e-6, max_weight=4)
        for b in backends(robs)
            @test termdict(propagate(rcirc, b, rthetas; min_abs_coeff=1e-6, max_weight=4)) ==
                  termdict(propagate(rcirc, b, rthetas; min_abs_coeff=1e-6, customtruncfunc=(mstr, coeff) -> get_weight(mstr) > 4))
        end
    end

    @testset "every backend agrees on wide Majorana strings" begin
        # 40 spinful sites are 160 Majorana modes, held in an NTupleInteger of three 64-bit words
        N = 40
        circ, thetas = hubbard_layer(bricklayertopology(N), N)
        circ, thetas = repeat(circ, 2), repeat(thetas, 2)
        obs = MajoranaSum(N, :nupndn, 20)
        @test MajoranaPropagation.termtype(obs) <: NTupleInteger

        results = propagateall(circ, obs, thetas; min_abs_coeff=1e-4, max_unpaired=2)
        @test agree(results)
        @test length(first(results)) > 100
    end

    @testset "thread=false == thread=true on a multi-task run" begin
        nx, ny = 5, 4
        N = nx * ny
        circ, thetas = hubbard_layer(rectangletopology(nx, ny), N)
        circ, thetas = repeat(circ, 4), repeat(thetas, 4)
        obs = MajoranaSum(N, :nupndn, 10)
        min_abs_coeff = 1e-6

        threaded = propagate(circ, VectorMajoranaSum(obs), thetas; min_abs_coeff, thread=true)
        serial = propagate(circ, VectorMajoranaSum(obs), thetas; min_abs_coeff, thread=false)
        zoned = propagate(circ, MultiMajoranaSum(VectorMajoranaSum(obs), 8), thetas; min_abs_coeff)
        dict = propagate(circ, obs, thetas; min_abs_coeff)

        @test length(threaded) > 16384  # this run spans several tasks
        @test termdict(threaded) == termdict(serial)
        @test termdict(zoned) == termdict(serial)
        @test termdict(dict) == termdict(serial)
    end

    @testset "path-property truncations" begin
        N = 6
        circ, thetas = hubbard_layer(bricklayertopology(N), N)
        circ, thetas = repeat(circ, 2), repeat(thetas, 2)
        obs = MajoranaSum(N, :nupndn, 3)

        # `max_freq` counts every cosine and sine application along a path, so it needs headroom
        for kw in ((; min_abs_coeff=0.0, max_freq=12), (; min_abs_coeff=0.0, max_sins=2), (; min_abs_coeff=1e-4, max_freq=30, max_sins=3))
            results = propagateall(circ, obs, thetas; kw...)
            @test agree(results)
            @test length(first(results)) > 1
        end

        # wrapped coefficients stay wrapped, counters included
        wrapped = [propagate(circ, wrapcoefficients(b, MajoranaFrequencyTracker), thetas; min_abs_coeff=1e-4, max_freq=12) for b in backends(obs)]
        @test all(MajoranaPropagation.coefftype(w) <: MajoranaFrequencyTracker for w in wrapped)
        @test agree(wrapped)
        @test length(first(wrapped)) > 1
    end

    @testset "imaginary-time rotations" begin
        nf = 5
        gam = dense_majoranas(nf)
        gates, betas = random_rotation_circuit(nf, 30; R=ImaginaryMajoranaRotation)
        betas .*= 0.2
        # the input already satisfies every truncation below, which the fused kernels take for granted
        rho = random_even_msum(nf, 12)
        truncate!(rho; min_abs_coeff=1e-5, max_weight=4)

        for max_weight in (Inf, 4), min_abs_coeff in (0.0, 1e-5)
            @test agree(propagateall(gates, rho, betas; min_abs_coeff, max_weight, heisenberg=false))
        end

        # against the dense oracle
        for b in backends(rho)
            evolved = propagate(gates[1:3], b, betas[1:3]; heisenberg=false, min_abs_coeff=0.0)
            K = I
            for (g, beta) in zip(gates[1:3], betas[1:3])
                K = exp(-beta / 2 * dense_string(g.ms_int, nf, gam)) * K
            end
            @test isapprox(dense_sum(MajoranaSum(evolved), gam), K * dense_sum(rho, gam) * K', atol=1e-10)
        end

        # fermionic imaginary-time gates, starting from the identity as a state does
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

        for normalize_coeffs in (true, false), truncate_each_mr in (nothing, true),
            kw in ((; min_abs_coeff=1e-4), (; min_abs_coeff=1e-5, max_unpaired=4), (; min_abs_coeff=0.0, max_weight=6))
            results = propagateall(circ, rho0, taus; kw..., heisenberg=false, normalize_coeffs, truncate_each_mr)
            @test agree(results)
            @test length(first(results)) > 1
        end
    end

    @testset "a large interacting circuit" begin
        # tails above 2 * _MIN_ELEMS_PER_TASK, so threaded runs merge with several tasks
        nx, ny = 5, 4
        nf = nx * ny
        circ = FermionicRotation[]
        thetas = Float64[]
        for pair in rectangletopology(nx, ny)
            push!(circ, FermionicRotation(:hop, pair))
            push!(thetas, 0.1)
            push!(circ, FermionicRotation(:nn, pair))
            push!(thetas, 0.1)
        end

        evolved = collect(backends(MajoranaSum(nf, :n, (nf + 1) ÷ 2)))
        for _ in 1:4
            for b in evolved
                propagate!(circ, b, thetas; min_abs_coeff=1e-6)
            end
            @test agree(evolved)
        end
        @test length(first(evolved)) > 100_000
    end

    @testset "propagate! into a pre-sized cache" begin
        N = 6
        circ, thetas = hubbard_layer(bricklayertopology(N), N)
        obs = MajoranaSum(N, :nupndn, 3)
        min_abs_coeff = 1e-4

        stock = VectorMajoranaSum(obs)
        caches = [PropagationCache(VectorMajoranaSum(obs)), PropagationCache(MultiMajoranaSum(VectorMajoranaSum(obs), 4))]
        foreach(cache -> resize!(cache, 10_000), caches)
        for _ in 1:3
            propagate!(circ, stock, thetas; min_abs_coeff)
            foreach(cache -> propagate!(circ, cache, thetas; min_abs_coeff), caches)
        end
        @test termdict(VectorMajoranaSum(caches[1])) == termdict(stock)
        @test termdict(MultiMajoranaSum(caches[2])) == termdict(stock)
    end

    @testset "MultiMajoranaSum" begin
        circ, thetas = hubbard_layer(bricklayertopology(4), 4)
        spinful = propagate(circ, MajoranaSum(4, :nupndn, 2), thetas; min_abs_coeff=0.0)
        spinless = random_even_msum(10, 50)
        fock = FockState(4, :checkerboard, true)

        for msum in (spinful, spinless), zones_of in (MajoranaSum, VectorMajoranaSum)
            m = MultiMajoranaSum(zones_of(msum), 4)
            @test nzones(m) == 4
            @test sum(zonesizes(m)) == length(m) == length(msum)
            @test MajoranaPropagation.is_spinful(m) == MajoranaPropagation.is_spinful(msum)
            @test MajoranaSum(m) == msum
            @test MajoranaSum(VectorMajoranaSum(m)) == msum
            @test all(getcoeff(m, ms) == coeff for (ms, coeff) in msum)

            # every string sits in the zone that owns it
            @test all(all(zoneof(m, ms) == zone_id for ms in terms(zone)) for (zone_id, zone) in enumerate(m.zones))

            # re-zoning keeps the terms
            @test MajoranaSum(MultiMajoranaSum(m, 8)) == msum
            @test nzones(MultiMajoranaSum(m, 8)) == 8

            @test MajoranaSum(unwrapcoefficients(wrapcoefficients(m, MajoranaFrequencyTracker))) == msum
            @test occursin("MultiMajoranaSum", sprint(show, m))
        end
        @test overlapwithfock(MultiMajoranaSum(VectorMajoranaSum(spinful), 4), fock) ≈ overlapwithfock(spinful, fock)

        # the zones are found by bit parities, which needs a power of two of them
        @test_throws ArgumentError MultiMajoranaSum(spinful, 3)
    end

    @testset "PropagationBase primitives on Majorana sums" begin
        msum = random_even_msum(8, 40)
        for b in backends(msum)
            @test termdict(mapcoeffs(c -> 2c, b)) == Dict(k => 2v for (k, v) in msum)
            @test termdict(filterterms(t -> get_weight(t) <= 2, b)) == Dict(k => v for (k, v) in msum if get_weight(k) <= 2)
            @test maxabscoeff(b) == maximum(abs, coefficients(msum))
        end

        vmsum = VectorMajoranaSum(msum)
        sortterms!(vmsum)
        @test issorted(terms(vmsum)) && termdict(vmsum) == termdict(msum)
    end

    # the results of the gates one by one and of the layers `tolayers` gathers them into, from every backend
    function flatandlayered(circ, thetas, obs; kwargs...)
        layered, layered_thetas = tolayers(circ, thetas)
        return [propagate(c, b, t; kwargs...) for b in backends(obs) for (c, t) in ((circ, thetas), (layered, layered_thetas))]
    end

    classplan(msum, layer, thetas) = MajoranaPropagation._classplan(layer, PropagationCache(msum), thetas, nothing)

    @testset "layers of commuting rotations propagate like their gates" begin
        N = 6
        circ, thetas = hubbard_layer(bricklayertopology(N), N)
        circ, thetas = repeat(circ, 2), repeat(thetas, 2)
        obs = MajoranaSum(N, :nupndn, 3)

        for min_abs_coeff in (0.0, 1e-4), max_weight in (Inf, 6), max_unpaired in (Inf, 2)
            results = flatandlayered(circ, thetas, obs; min_abs_coeff, max_weight, max_unpaired)
            @test agree(results)
            @test length(first(results)) > 1
        end

        # a hopping truncates only after its last rotation, or after each, and its layer does the same
        default, each = (flatandlayered(circ, thetas, obs; min_abs_coeff=1e-5, max_unpaired=2, truncate_each_mr) for truncate_each_mr in (nothing, true))
        @test agree(default) && agree(each)
        @test termdict(first(default)) != termdict(first(each))

        # the layers are applied class by class on every backend, unless the coefficients track frequencies
        layered, layered_thetas = tolayers(circ, thetas)
        layers = [(layer, theta) for (layer, theta) in zip(layered, layered_thetas) if layer isa GateLayer]
        @test length(layers) == 6
        for b in backends(obs)
            @test all(!isnothing(classplan(b, layer, theta)) for (layer, theta) in layers)
            @test all(isnothing(classplan(wrapcoefficients(b, MajoranaFrequencyTracker), layer, theta)) for (layer, theta) in layers)
        end
        @test agree(flatandlayered(circ, thetas, obs; min_abs_coeff=1e-4, max_freq=12))

        # random Majorana rotations that commute, some of them on the same Majoranas
        nf = 8
        rotations, angles = random_rotation_circuit(nf, 60)
        kept = Int[]
        for i in eachindex(rotations)
            if all(j -> MajoranaPropagation.commutes(rotations[i].ms_int, rotations[j].ms_int), kept)
                push!(kept, i)
            end
        end
        layer = GateLayer(rotations[kept]; guaranteed_commutes=true)
        robs = random_even_msum(nf, 5)
        truncate!(robs; min_abs_coeff=1e-6, max_weight=4)
        @test length(kept) > 4 && !isnothing(classplan(robs, layer, angles[kept]))
        for max_weight in (Inf, 4), min_abs_coeff in (0.0, 1e-6)
            @test agree([propagate(c, b, t; min_abs_coeff, max_weight) for b in backends(robs) for (c, t) in ((rotations[kept], angles[kept]), (layer, angles[kept]))])
        end

        # rotations that do not commute, forced into a layer, apply one by one
        forced = GateLayer(circ[1:8]; guaranteed_commutes=true)
        @test isnothing(classplan(obs, forced, thetas[1:8]))
        @test termdict(propagate(forced, obs, thetas[1:8]; min_abs_coeff=1e-4)) == termdict(propagate(circ[1:8], obs, thetas[1:8]; min_abs_coeff=1e-4))

        # wide strings, held in NTupleIntegers of three words
        circ40, thetas40 = hubbard_layer(bricklayertopology(40), 40)
        @test agree(flatandlayered(repeat(circ40, 2), repeat(thetas40, 2), MajoranaSum(40, :nupndn, 20); min_abs_coeff=1e-4, max_unpaired=2))

        # a vector sum large enough for several tasks
        circ54, thetas54 = hubbard_layer(rectangletopology(5, 4), 20)
        circ54, thetas54 = repeat(circ54, 4), repeat(thetas54, 4)
        layered54, layered_thetas54 = tolayers(circ54, thetas54)
        obs54 = VectorMajoranaSum(MajoranaSum(20, :nupndn, 10))
        flat54 = propagate(circ54, obs54, thetas54; min_abs_coeff=1e-6)
        @test length(flat54) > 16384
        @test termdict(propagate(layered54, obs54, layered_thetas54; min_abs_coeff=1e-6)) == termdict(flat54)

        # imaginary-time rotations normalize the whole sum, so their layers apply them one by one
        icirc = [ImaginaryFermionicRotation(:nupndn, [i]) for i in 1:N]
        ilayered, itaus = tolayers(icirc, fill(0.2, N))
        rho0 = MajoranaSum(Float64, N, true)
        add!(rho0, getinttype(2N)(0), 1.0)
        @test only(ilayered) isa GateLayer
        @test termdict(propagate(ilayered, rho0, itaus; heisenberg=false)) == termdict(propagate(icirc, rho0, fill(0.2, N); heisenberg=false))
    end

    @testset "tolayers gathers the rotations that commute" begin
        qinds = PauliPropagation.PropagationBase.qinds
        @test qinds(MajoranaRotation(MajoranaString(4, [1, 2, 5, 8]))) == [1, 2, 5, 8]
        # the Majoranas of the up and the down hopping of a bond differ
        @test isdisjoint(qinds(FermionicRotation(:hopup, [1, 2])), qinds(FermionicRotation(:hopdn, [1, 2])))
        @test PauliPropagation.commutes(FermionicRotation(:nn, [1, 2]), FermionicRotation(:nn, [2, 3]))
        @test !PauliPropagation.commutes(FermionicRotation(:hop, [1, 2]), FermionicRotation(:hop, [2, 3]))
        @test length(tolayers([FermionicRotation(:nn, pair) for pair in staircasetopology(6)])) == 1
        @test length(tolayers([FermionicRotation(:hop, pair) for pair in staircasetopology(6)])) == 5
    end
end
