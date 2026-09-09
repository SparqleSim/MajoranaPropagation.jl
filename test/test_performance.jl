# Fused truncations (src/Performance): every overload must reproduce the default path exactly.
using MajoranaPropagation
using MajoranaPropagation.Performance
using PauliPropagation
using Test
using Random

if !isdefined(Main, :dense_majoranas)
    include("testutils_dense.jl")
end

const perf_rng = MersenneTwister(7)

# one first-order Trotter layer of the spinful Hubbard model: hoppings decompose into two rotations
# that are only number-preserving together, `:nupndn` into rotations truncated one by one
function perf_hubbard_layer(topo, n_sites; t=1.0, U=1.0, dt=0.1)
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

# random even-weight Majorana rotations on `nf` spinless fermions
function perf_rotation_circuit(nf, n_gates)
    circ = MajoranaRotation{getinttype(nf)}[]
    thetas = Float64[]
    for _ in 1:n_gates
        w = rand(perf_rng, (2, 4))
        gammas = Int[]
        while length(gammas) < w
            g = rand(perf_rng, 1:2*nf)
            g in gammas || push!(gammas, g)
        end
        push!(circ, MajoranaRotation(MajoranaString(nf, gammas)))
        push!(thetas, randn(perf_rng))
    end
    return circ, thetas
end

# O(N) exact comparison of two sums of any backend
perf_dict(msum) = Dict(k => v for (k, v) in msum)

@testset "Performance" begin

    @testset "opt-in: fused=false is the default path" begin
        N = 6
        circ, thetas = perf_hubbard_layer(bricklayertopology(N), N)
        circ, thetas = repeat(circ, 2), repeat(thetas, 2)
        obs = MajoranaSum(N, :nupndn, 3)
        min_abs_coeff = 1e-4

        stock_dict = propagate(circ, obs, thetas; min_abs_coeff)
        stock_vec = propagate(circ, VectorMajoranaSum(obs), thetas; min_abs_coeff)

        @test Performance.propagate(circ, obs, thetas; min_abs_coeff, fused=false) == stock_dict
        @test perf_dict(Performance.propagate(circ, VectorMajoranaSum(obs), thetas; min_abs_coeff, fused=false)) == perf_dict(stock_vec)

        # a fused call must not change what a later plain call returns
        Performance.propagate(circ, obs, thetas; min_abs_coeff)
        Performance.propagate(circ, VectorMajoranaSum(obs), thetas; min_abs_coeff)
        @test propagate(circ, obs, thetas; min_abs_coeff) == stock_dict
        @test perf_dict(propagate(circ, VectorMajoranaSum(obs), thetas; min_abs_coeff)) == perf_dict(stock_vec)
    end

    @testset "fused vector == default under every truncation" begin
        N = 6
        circ, thetas = perf_hubbard_layer(bricklayertopology(N), N)
        circ, thetas = repeat(circ, 2), repeat(thetas, 2)
        obs = MajoranaSum(N, :nupndn, 3)

        for min_abs_coeff in (0.0, 1e-4), max_weight in (Inf, 4, 6), max_unpaired in (Inf, 2, 4)
            kw = (; min_abs_coeff, max_weight, max_unpaired)
            stock = propagate(circ, VectorMajoranaSum(obs), thetas; kw...)
            fused = Performance.propagate(circ, VectorMajoranaSum(obs), thetas; kw...)
            @test perf_dict(fused) == perf_dict(stock)
            @test length(fused) > 1
        end

        # both ways of truncating a fermionic gate's rotations
        for truncate_each_mr in (true, false)
            kw = (; min_abs_coeff=1e-4, max_unpaired=4, truncate_each_mr)
            stock = propagate(circ, VectorMajoranaSum(obs), thetas; kw...)
            fused = Performance.propagate(circ, VectorMajoranaSum(obs), thetas; kw...)
            @test perf_dict(fused) == perf_dict(stock)
        end

        # a precomputed unpaired mask and a custom truncation
        let mask = create_unpaired_mask(2N), customtruncfunc = (mstr, coeff) -> get_weight(mstr) == 4 && abs(coeff) < 1e-2
            kw = (; min_abs_coeff=1e-5, max_unpaired=4, unpaired_mask=mask, customtruncfunc)
            stock = propagate(circ, VectorMajoranaSum(obs), thetas; kw...)
            fused = Performance.propagate(circ, VectorMajoranaSum(obs), thetas; kw...)
            @test perf_dict(fused) == perf_dict(stock)
        end
    end

    @testset "fused vector == default for bare Majorana rotations" begin
        nf = 8
        circ, thetas = perf_rotation_circuit(nf, 40)
        obs = random_even_msum(nf, 5)
        truncate!(obs; min_abs_coeff=1e-6, max_weight=4)

        for max_weight in (Inf, 4), min_abs_coeff in (0.0, 1e-6)
            kw = (; min_abs_coeff, max_weight)
            stock = propagate(circ, VectorMajoranaSum(obs), thetas; kw...)
            fused = Performance.propagate(circ, VectorMajoranaSum(obs), thetas; kw...)
            @test perf_dict(fused) == perf_dict(stock)
        end
    end

    @testset "fused vector == default on wide Majorana strings" begin
        # 40 spinful sites are 160 Majorana modes, well beyond a machine word
        N = 40
        circ, thetas = perf_hubbard_layer(bricklayertopology(N), N)
        circ, thetas = repeat(circ, 2), repeat(thetas, 2)
        obs = MajoranaSum(N, :nupndn, 20)
        @test sizeof(MajoranaPropagation.majoranatype(obs)) > 8

        kw = (; min_abs_coeff=1e-4, max_unpaired=2)
        stock = propagate(circ, VectorMajoranaSum(obs), thetas; kw...)
        fused = Performance.propagate(circ, VectorMajoranaSum(obs), thetas; kw...)
        @test perf_dict(fused) == perf_dict(stock)
        @test length(fused) > 100
    end

    @testset "fused vector: thread=false == thread=true on a multi-task run" begin
        nx, ny = 5, 4
        N = nx * ny
        circ, thetas = perf_hubbard_layer(rectangletopology(nx, ny), N; dt=0.1)
        circ, thetas = repeat(circ, 4), repeat(thetas, 4)
        obs = MajoranaSum(N, :nupndn, 10)
        min_abs_coeff = 1e-6

        threaded = Performance.propagate(circ, VectorMajoranaSum(obs), thetas; min_abs_coeff, thread=true)
        serial = Performance.propagate(circ, VectorMajoranaSum(obs), thetas; min_abs_coeff, thread=false)
        stock = propagate(circ, VectorMajoranaSum(obs), thetas; min_abs_coeff)

        @test length(threaded) > 16384  # sanity check that this run actually spans several tasks
        @test perf_dict(threaded) == perf_dict(serial)
        @test perf_dict(threaded) == perf_dict(stock)
    end

    @testset "fused vector == default for path-property truncations" begin
        N = 6
        circ, thetas = perf_hubbard_layer(bricklayertopology(N), N)
        circ, thetas = repeat(circ, 2), repeat(thetas, 2)
        obs = MajoranaSum(N, :nupndn, 3)

        # `max_freq` counts every cosine and sine application along a path, so it needs headroom
        for kw in ((; min_abs_coeff=0.0, max_freq=12), (; min_abs_coeff=0.0, max_sins=2), (; min_abs_coeff=1e-4, max_freq=30, max_sins=3))
            stock = propagate(circ, VectorMajoranaSum(obs), thetas; kw...)
            fused = Performance.propagate(circ, VectorMajoranaSum(obs), thetas; kw...)
            @test perf_dict(fused) == perf_dict(stock)
            @test length(fused) > 1
        end

        # wrapped coefficients stay wrapped, counters included
        wrapped = wrapcoefficients(VectorMajoranaSum(obs), MajoranaFrequencyTracker)
        stock = propagate(circ, wrapped, thetas; min_abs_coeff=1e-4, max_freq=12)
        fused = Performance.propagate(circ, wrapped, thetas; min_abs_coeff=1e-4, max_freq=12)
        @test MajoranaPropagation.coefftype(fused) <: MajoranaFrequencyTracker
        @test length(fused) > 1
        @test perf_dict(fused) == perf_dict(stock)
    end

    @testset "fused vector == default for imaginary-time rotations" begin
        nf = 5
        gam = dense_majoranas(nf)
        gates = ImaginaryMajoranaRotation{getinttype(nf)}[]
        betas = Float64[]
        for _ in 1:30
            push!(gates, ImaginaryMajoranaRotation(MajoranaString(nf, Int(random_even_string(nf)))))
            push!(betas, 0.2 * randn(perf_rng))
        end
        rho = random_even_msum(nf, 4)

        for max_weight in (Inf, 4), min_abs_coeff in (0.0, 1e-5)
            kw = (; min_abs_coeff, max_weight, heisenberg=false)
            stock = propagate(gates, VectorMajoranaSum(rho), betas; kw...)
            fused = Performance.propagate(gates, VectorMajoranaSum(rho), betas; kw...)
            @test perf_dict(fused) == perf_dict(stock)
        end

        # against the dense oracle, as for the default path
        evolved = Performance.propagate(gates[1:3], VectorMajoranaSum(rho), betas[1:3]; heisenberg=false, min_abs_coeff=0.0)
        K = I
        for (g, b) in zip(gates[1:3], betas[1:3])
            K = exp(-b / 2 * dense_string(g.ms_int, nf, gam)) * K
        end
        @test isapprox(dense_sum(MajoranaSum(nf, false, perf_dict(evolved)), gam), K * dense_sum(rho, gam) * K', atol=1e-10)

        # fermionic imaginary-time gates, starting from the identity as a state does
        N = 4
        topo = bricklayertopology(N)
        circ = ImaginaryFermionicRotation[]
        taus = Float64[]
        for pair in topo
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

        for normalize_coeffs in (true, false), truncate_each_mr in (nothing, true), kw in ((; min_abs_coeff=1e-4), (; min_abs_coeff=1e-5, max_unpaired=4), (; min_abs_coeff=0.0, max_weight=6))
            kw = (; kw..., heisenberg=false, normalize_coeffs, truncate_each_mr)
            stock = propagate(circ, VectorMajoranaSum(rho0), taus; kw...)
            fused = Performance.propagate(circ, VectorMajoranaSum(rho0), taus; kw...)
            @test perf_dict(fused) == perf_dict(stock)
            @test length(fused) > 1
        end
    end

    @testset "fused Dict == default under every truncation" begin
        N = 6
        circ, thetas = perf_hubbard_layer(bricklayertopology(N), N)
        circ, thetas = repeat(circ, 2), repeat(thetas, 2)
        obs = MajoranaSum(N, :nupndn, 3)

        for min_abs_coeff in (0.0, 1e-3, 1e-5), max_weight in (Inf, 4, 6), max_unpaired in (Inf, 2, 4), truncate_each_mr in (nothing, false)
            kw = (; min_abs_coeff, max_weight, max_unpaired, truncate_each_mr)
            stock = propagate(circ, obs, thetas; kw...)
            fused = Performance.propagate(circ, obs, thetas; kw...)
            @test fused == stock
        end

        # bare rotations too
        nf = 8
        rcirc, rthetas = perf_rotation_circuit(nf, 40)
        robs = random_even_msum(nf, 5)
        truncate!(robs; min_abs_coeff=1e-6, max_weight=4)
        for kw in ((; min_abs_coeff=0.0, max_weight=4), (; min_abs_coeff=1e-6), (; min_abs_coeff=1e-4, max_weight=4))
            @test Performance.propagate(rcirc, robs, rthetas; kw...) == propagate(rcirc, robs, rthetas; kw...)
        end

        # the Dict and vector fused paths agree with each other as well
        stock = propagate(circ, obs, thetas; min_abs_coeff=1e-4, max_unpaired=4)
        fused_vec = Performance.propagate(circ, VectorMajoranaSum(obs), thetas; min_abs_coeff=1e-4, max_unpaired=4)
        @test perf_dict(fused_vec) == perf_dict(stock)
    end

    @testset "propagate! into a cache" begin
        N = 6
        circ, thetas = perf_hubbard_layer(bricklayertopology(N), N)
        obs = MajoranaSum(N, :nupndn, 3)
        min_abs_coeff = 1e-4

        cache = PropagationCache(VectorMajoranaSum(obs))
        resize!(cache, 10_000)
        stock = VectorMajoranaSum(obs)
        for _ in 1:3
            Performance.propagate!(circ, cache, thetas; min_abs_coeff)
            propagate!(circ, stock, thetas; min_abs_coeff)
        end
        @test perf_dict(VectorMajoranaSum(cache)) == perf_dict(stock)
    end
end
