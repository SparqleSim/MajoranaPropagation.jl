# rewindgradient: the expectation value is that of `propagate`, and the gradient agrees with central finite
# differences of it, for bare Majorana rotations, fermionic rotations and layers of them, from every backend.
# Without truncation the backward sweep undoes the forward sweep exactly, so the comparisons run at
# `min_abs_coeff=0.0`, which truncates nothing.
using MajoranaPropagation
using PauliPropagation
using Test
using Random

const gradient_rng = MersenneTwister(11)

# the derivative of overlapfunc(propagate(circ, obs, thetas)) by every angle, by central finite differences
function finitedifferences(circ, obs, thetas, overlapfunc; step=1e-6, kwargs...)
    shifted(k, by) = [i == k ? theta + by : theta for (i, theta) in enumerate(thetas)]
    expec(angles) = overlapfunc(propagate(circ, obs, angles; kwargs...))
    return [(expec(shifted(k, step)) - expec(shifted(k, -step))) / 2step for k in eachindex(thetas)]
end

function agreeswithfinitedifferences(circ, obs, thetas, overlapfunc)
    expec, grad = rewindgradient(circ, obs, thetas, overlapfunc; min_abs_coeff=0.0)
    return expec ≈ overlapfunc(propagate(circ, obs, thetas; min_abs_coeff=0.0)) &&
           isapprox(grad, finitedifferences(circ, obs, thetas, overlapfunc; min_abs_coeff=0.0); atol=1e-7)
end

gradient_backends(obs) = (obs, VectorMajoranaSum(obs), MultiMajoranaSum(obs, 4), MultiMajoranaSum(VectorMajoranaSum(obs), 4))

@testset "rewindgradient" begin
    # spinless: hoppings, interactions and pairings, which are not number-preserving, a bare Majorana rotation
    # that does not commute with its neighbours, and an observable of several strings
    n = 4
    obs = MajoranaSum(n, :n, 2)
    add!(obs, MajoranaSum(n, :hop, (1, 3)))
    TT = getinttype(n)
    fock = FockState(n, [1, 3])
    overlapfunc(msum) = overlapwithfock(msum, fock)
    circ = PauliPropagation.Gate[]
    for site in 1:3
        for pair in ((1, 2), (3, 4), (2, 3))
            append!(circ, [FermionicRotation(:hop, pair), FermionicRotation(:nn, pair), FermionicRotation(:pair, pair)])
        end
        append!(circ, [MajoranaRotation(TT(0b00111100)), FermionicRotation(:n, site)])
    end
    thetas = randn(gradient_rng, length(circ))
    expec, grad = rewindgradient(circ, obs, thetas, overlapfunc; min_abs_coeff=0.0)

    @testset "bare Majorana rotations" begin
        rotations = [MajoranaRotation(TT(mask)) for mask in (0b00000011, 0b00000110, 0b00111100, 0b11000011, 0b01100110, 0b00001111)]
        @test agreeswithfinitedifferences(rotations, obs, randn(gradient_rng, length(rotations)), overlapfunc)
        # a single gate with its angle is a circuit of one gate
        @test rewindgradient(rotations[3], obs, 0.4, overlapfunc; min_abs_coeff=0.0) ==
              rewindgradient([rotations[3]], obs, [0.4], overlapfunc; min_abs_coeff=0.0)
    end

    @testset "fermionic rotations" begin
        @test agreeswithfinitedifferences(circ, obs, thetas, overlapfunc)
    end

    @testset "a spinful Hubbard layer" begin
        n_sites = 3
        hubbard = FermionicRotation[]
        for pair in bricklayertopology(n_sites)
            append!(hubbard, [FermionicRotation(:hopup, pair), FermionicRotation(:hopdn, pair)])
        end
        append!(hubbard, [FermionicRotation(:nupndn, site) for site in 1:n_sites])
        doublons = MajoranaSum(n_sites, :nupndn, 2)
        checkerboard = FockState(n_sites, :checkerboard, true)
        @test agreeswithfinitedifferences([hubbard; hubbard], doublons, randn(gradient_rng, 2 * length(hubbard)),
            msum -> overlapwithfock(msum, checkerboard))
    end

    @testset "every backend, and a cache" begin
        for msum in gradient_backends(obs)
            backend_expec, backend_grad = rewindgradient(circ, msum, thetas, overlapfunc; min_abs_coeff=0.0)
            @test backend_expec ≈ expec
            @test backend_grad ≈ grad
        end

        # the in-place version leaves the propagated operator in the cache
        cache = PropagationCache(VectorMajoranaSum(obs))
        cache_expec, cache_grad = rewindgradient!(circ, cache, thetas, overlapfunc; min_abs_coeff=0.0)
        @test cache_expec ≈ expec
        @test cache_grad ≈ grad
        @test MajoranaSum(cache) == MajoranaSum(propagate(circ, VectorMajoranaSum(obs), thetas; min_abs_coeff=0.0))
    end

    @testset "layers" begin
        layered, layered_thetas = torotationlayers(circ, thetas)
        @test any(gate -> gate isa GateLayer, layered)
        for msum in gradient_backends(obs)
            layered_expec, layered_grad = rewindgradient(layered, msum, layered_thetas, overlapfunc; min_abs_coeff=0.0)
            @test layered_expec ≈ expec
            # a vector of derivatives for every layer, as its angles are a vector
            @test all(length(layer_grad) == length(layer_thetas) for (layer_grad, layer_thetas) in zip(layered_grad, layered_thetas))
            @test reduce(vcat, layered_grad) ≈ grad
        end
    end

    @testset "overlap with a Majorana sum" begin
        target = propagate(circ[1:6], MajoranaSum(n, :nn, (2, 3)), randn(gradient_rng, 6); min_abs_coeff=0.0)
        @test agreeswithfinitedifferences(circ, obs, thetas, msum -> MajoranaPropagation.scalarproduct(msum, target))
    end

    @testset "truncations" begin
        # nothing is truncated, so it does not matter which strings the dual keeps
        @test last(rewindgradient(circ, obs, thetas, overlapfunc; min_abs_coeff=0.0, keep_rel_dual=Inf)) ≈ grad
        # the default threshold drops what the backward sweep cancels, and nothing else
        @test last(rewindgradient(circ, obs, thetas, overlapfunc)) ≈ grad

        # a custom truncation reads the coefficient of the operator in both sweeps
        keepsall(mstr, coeff) = !(coeff isa Float64)
        @test last(rewindgradient(circ, obs, thetas, overlapfunc; min_abs_coeff=0.0, customtruncfunc=keepsall)) ≈ grad

        for msum in gradient_backends(obs), kwargs in ((; min_abs_coeff=1e-2), (; min_abs_coeff=1e-2, keep_rel_dual=Inf),
            (; min_rel_coeff=1e-2), (; max_weight=4), (; max_unpaired=2))
            truncated_expec, truncated_grad = rewindgradient(circ, msum, thetas, overlapfunc; kwargs...)
            @test truncated_expec ≈ overlapfunc(propagate(circ, msum, thetas; filter(pair -> first(pair) != :keep_rel_dual, pairs(kwargs))...))
            @test length(truncated_grad) == length(thetas) && all(isfinite, truncated_grad)
        end

        # truncating everything leaves no expectation and no gradient
        empty_expec, empty_grad = rewindgradient(circ, obs, thetas, overlapfunc; min_abs_coeff=10.0)
        @test empty_expec == 0 && all(iszero, empty_grad)
    end

    @testset "what cannot be undone is rejected" begin
        imaginary = [ImaginaryMajoranaRotation(TT(0b00000011))]
        @test_throws ArgumentError rewindgradient(imaginary, obs, [0.3], overlapfunc)
        frozen = [circ[1], freeze(circ[2], 0.3)]
        @test_throws ArgumentError rewindgradient(frozen, obs, [0.3], overlapfunc)
        tracked = wrapcoefficients(obs, MajoranaFrequencyTracker)
        @test_throws ArgumentError rewindgradient(circ, tracked, thetas, overlapfunc)
    end
end
