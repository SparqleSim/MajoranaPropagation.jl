using MajoranaPropagation
using MajoranaPropagation.GateLookup
using PauliPropagation   # for PauliRotation, used by the Yao helper
using Test
using Random

Random.seed!(1234)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# coefficient dict keyed by BigInt so sums with different integer widths compare cleanly
function _coeffdict(msum)
    d = Dict{BigInt,ComplexF64}()
    for (k, v) in msum
        d[BigInt(k)] = ComplexF64(v)
    end
    return d
end

function _msums_match(a, b; atol=1e-9)
    da, db = _coeffdict(a), _coeffdict(b)
    for k in union(keys(da), keys(db))
        if abs(get(da, k, 0.0im) - get(db, k, 0.0im)) > atol
            return false
        end
    end
    return true
end

# propagate `obs` through the lookup gate and through the reference FermionicRotation, compare
function _lookup_matches_fr(symbol, is_spinful, site_inds, theta, obs; atol=1e-9)
    fr = FermionicRotation(symbol, site_inds)
    ref = propagate([fr], deepcopy(obs), [theta]; min_abs_coeff=1e-14)

    look = FermionicRotationLookup(symbol, length(site_inds), is_spinful, site_inds, theta)
    got = propagate([look], deepcopy(obs); min_abs_coeff=1e-14)

    return _msums_match(ref, got; atol=atol)
end

# build a FermionicRotation circuit and the matching FermionicRotationLookup circuit from a list
# of `(symbol, sites, theta)` specs, propagate `obs` through both, and compare. This exercises
# applying many gates one after the other.
function _circuit_matches_fr(specs, is_spinful, obs; atol=1e-8, min_abs_coeff=1e-12)
    fr_circ = FermionicRotation[]
    look_circ = []
    thetas = Float64[]
    for (sym, sites, th) in specs
        push!(fr_circ, FermionicRotation(sym, sites))
        push!(look_circ, FermionicRotationLookup(sym, length(sites), is_spinful, sites, th))
        push!(thetas, th)
    end
    ref = propagate(fr_circ, deepcopy(obs), thetas; min_abs_coeff=min_abs_coeff)
    got = propagate(look_circ, deepcopy(obs); min_abs_coeff=min_abs_coeff)
    return _msums_match(ref, got; atol=atol)
end

# ---------------------------------------------------------------------------
# MajoranaTransferMap structure
# ---------------------------------------------------------------------------

@testset "MajoranaTransferMap structure" begin
    # spinless 2-site gate => 4 Majorana modes => 2^4 = 16 columns
    gate = FermionicRotationLookup(:hop, 2, false, [1, 2], 0.4)
    tmap = gate.transfer_map
    @test MajoranaPropagation.GateLookup.ncolumns(tmap) == 2^4
    # the identity column (input string 0) maps to itself: cumulative string 0, coefficient 1
    @test length(tmap[0]) == 1
    @test tmap[0][1] == (0x00, 1.0 + 0.0im)
    @test_throws BoundsError tmap[-1]
    @test_throws BoundsError tmap[2^4]

    # spinful 1-site gate => 4 Majorana modes => 16 columns
    sgate = FermionicRotationLookup(:nupndn, 1, true, [1], 0.4)
    @test MajoranaPropagation.GateLookup.ncolumns(sgate.transfer_map) == 2^4

    # spinful 2-site gate => 8 Majorana modes => 256 columns
    sgate2 = FermionicRotationLookup(:hopup, 2, true, [1, 2], 0.4)
    @test MajoranaPropagation.GateLookup.ncolumns(sgate2.transfer_map) == 2^8
end

# ---------------------------------------------------------------------------
# spinless: lookup gate matches FermionicRotation
# ---------------------------------------------------------------------------

@testset "spinless lookup == FermionicRotation" begin
    n_sites = 6
    thetas = [0.0, 0.13, 0.4, 1.0, -0.7, 2.4]

    # observables, including ones that span the gate boundary (odd restricted weight)
    observables = [
        MajoranaSum(n_sites, :n, 3),
        MajoranaSum(n_sites, :nn, [3, 4]),
        MajoranaSum(n_sites, :hop, [3, 5]),            # spans gate boundary
        MajoranaSum(n_sites, :hop, [1, 6]),
        MajoranaSum(n_sites, :pair, [2, 4]),
        MajoranaSum(n_sites, :n, 2) + MajoranaSum(n_sites, :hop, [2, 5]),
    ]

    gateconfigs = [
        (:n, [3]),
        (:hop, [2, 3]),
        (:hop, [3, 5]),    # non-contiguous
        (:nn, [3, 4]),
        (:nn, [2, 5]),     # non-contiguous
        (:pair, [3, 4]),
    ]

    @testset "symbol=$symbol sites=$sites" for (symbol, sites) in gateconfigs
        for theta in thetas, obs in observables
            @test _lookup_matches_fr(symbol, false, sites, theta, obs)
        end
    end
end

# ---------------------------------------------------------------------------
# spinful: lookup gate matches FermionicRotation
# ---------------------------------------------------------------------------

@testset "spinful lookup == FermionicRotation" begin
    n_sites = 4
    thetas = [0.0, 0.21, 0.85, -0.55, 1.7]

    observables = [
        MajoranaSum(n_sites, :nup, 2),
        MajoranaSum(n_sites, :ndn, 3),
        MajoranaSum(n_sites, :nupndn, 2),
        MajoranaSum(n_sites, :hole, 3),
        MajoranaSum(n_sites, :hopup, [2, 4]),          # spans gate boundary
        MajoranaSum(n_sites, :hopdn, [1, 3]),
        MajoranaSum(n_sites, :nup, 2) + MajoranaSum(n_sites, :ndn, 3),
    ]

    gateconfigs = [
        (:nup, [2]),
        (:ndn, [3]),
        (:nupndn, [2]),
        (:hole, [3]),
        (:hopup, [2, 3]),
        (:hopdn, [2, 4]),    # non-contiguous
        (:pairup, [2, 3]),
        (:pairdn, [1, 3]),   # non-contiguous
    ]

    @testset "symbol=$symbol sites=$sites" for (symbol, sites) in gateconfigs
        for theta in thetas, obs in observables
            @test _lookup_matches_fr(symbol, true, sites, theta, obs)
        end
    end
end

# ---------------------------------------------------------------------------
# arbitrary site ordering and direction-sensitive gates
# ---------------------------------------------------------------------------

@testset "arbitrary site order" begin
    n_sites = 6
    # symmetric gates with site indices passed in descending / arbitrary order must agree with
    # FermionicRotation (which sorts internally)
    observables_sl = [MajoranaSum(n_sites, :n, 3), MajoranaSum(n_sites, :hop, [3, 5])]
    for obs in observables_sl, theta in [0.3, -0.9]
        @test _lookup_matches_fr(:hop, false, [2, 3], theta, obs)
        @test _lookup_matches_fr(:nn, false, [2, 5], theta, obs)
        @test _lookup_matches_fr(:pair, false, [1, 4], theta, obs)
    end

    # direction-sensitive spinful gate: :hopupdn(up@s1, down@s2) is NOT symmetric in its sites,
    # so both orderings must match their respective FermionicRotation
    observables_sf = [MajoranaSum(n_sites, :nup, 2), MajoranaSum(n_sites, :nupndn, 4),
        MajoranaSum(n_sites, :hopup, [2, 5])]
    for obs in observables_sf, theta in [0.4, -1.1]
        @test _lookup_matches_fr(:hopupdn, true, [2, 4], theta, obs)
        @test _lookup_matches_fr(:hopupdn, true, [4, 2], theta, obs)   # reversed up/down sites
        @test _lookup_matches_fr(:hopupdn, true, [5, 1], theta, obs)
    end
end

# ---------------------------------------------------------------------------
# repeated site indices
# ---------------------------------------------------------------------------

@testset "repeated site indices" begin
    # `site_inds` may repeat a site whenever FermionicRotation(symbol, site_inds) is a valid
    # operator with that repetition. The lookup gate must build without error and still match.
    n_sites = 5
    observables_sl = [MajoranaSum(n_sites, :n, 3), MajoranaSum(n_sites, :hop, [2, 4])]
    for obs in observables_sl, theta in [0.4, -0.8]
        @test _lookup_matches_fr(:nn, false, [2, 2], theta, obs)
        @test _lookup_matches_fr(:hop, false, [3, 3], theta, obs)
        @test _lookup_matches_fr(:pair, false, [2, 2], theta, obs)
    end

    observables_sf = [MajoranaSum(n_sites, :nup, 3), MajoranaSum(n_sites, :nupndn, 2)]
    for obs in observables_sf, theta in [0.4, -0.8]
        @test _lookup_matches_fr(:hopup, true, [2, 2], theta, obs)
        @test _lookup_matches_fr(:hopdn, true, [3, 3], theta, obs)
        @test _lookup_matches_fr(:pairup, true, [2, 2], theta, obs)
    end
end

# ---------------------------------------------------------------------------
# complex-coefficient observables
# ---------------------------------------------------------------------------

@testset "complex observables" begin
    n_sites = 4
    for obs in [MajoranaSum(n_sites, :f, 2), MajoranaSum(n_sites, :fdag, 3)]
        @test _lookup_matches_fr(:hop, false, [1, 2], 0.33, obs)
        @test _lookup_matches_fr(:nn, false, [2, 3], 0.77, obs)
    end
end

# ---------------------------------------------------------------------------
# multi-gate circuits
# ---------------------------------------------------------------------------

@testset "multi-gate circuit (spinless)" begin
    n_sites = 6
    topo = bricklayertopology(n_sites)
    obs = MajoranaSum(n_sites, :n, 3) + MajoranaSum(n_sites, :nn, [2, 5])

    fr_circ = FermionicRotation[]
    look_circ = []
    thetas = Float64[]
    for (i, j) in topo
        push!(fr_circ, FermionicRotation(:hop, [i, j]))
        push!(look_circ, FermionicRotationLookup(:hop, 2, false, [i, j], 0.2))
        push!(thetas, 0.2)
    end
    for (i, j) in topo
        push!(fr_circ, FermionicRotation(:nn, [i, j]))
        push!(look_circ, FermionicRotationLookup(:nn, 2, false, [i, j], 0.5))
        push!(thetas, 0.5)
    end

    ref = propagate(fr_circ, deepcopy(obs), thetas; min_abs_coeff=1e-12)
    got = propagate(look_circ, deepcopy(obs); min_abs_coeff=1e-12)
    @test length(ref) == length(got)
    @test _msums_match(ref, got; atol=1e-8)
end

@testset "multi-gate circuit (spinful)" begin
    n_sites = 4
    obs = MajoranaSum(n_sites, :nupndn, 2)

    fr_circ = FermionicRotation[]
    look_circ = []
    thetas = Float64[]
    for (sym, sites, th) in [(:hopup, [1, 2], 0.3), (:hopdn, [2, 3], 0.25),
        (:hopup, [3, 4], 0.4), (:nupndn, [2], 0.6), (:nupndn, [3], 0.6)]
        push!(fr_circ, FermionicRotation(sym, sites))
        push!(look_circ, FermionicRotationLookup(sym, length(sites), true, sites, th))
        push!(thetas, th)
    end

    ref = propagate(fr_circ, deepcopy(obs), thetas; min_abs_coeff=1e-12)
    got = propagate(look_circ, deepcopy(obs); min_abs_coeff=1e-12)
    @test _msums_match(ref, got; atol=1e-8)
end

# ---------------------------------------------------------------------------
# sequential multi-gate circuits (many gates one after the other)
# ---------------------------------------------------------------------------

@testset "sequential gates (spinless, mixed)" begin
    n_sites = 7
    # a hand-built circuit mixing gate types, non-contiguous sites, descending order, repeats
    specs = [
        (:hop, [1, 2], 0.31), (:hop, [3, 2], -0.4), (:nn, [2, 5], 0.7),
        (:pair, [4, 6], 0.2), (:hop, [6, 7], -0.55), (:nn, [1, 3], 0.9),
        (:hop, [2, 5], 0.15), (:nn, [4, 4], 0.5), (:pair, [1, 7], -0.25),
        (:hop, [5, 4], 0.6), (:nn, [3, 6], -0.8), (:hop, [7, 1], 0.33),
    ]
    for obs in [MajoranaSum(n_sites, :n, 4),
        MajoranaSum(n_sites, :hop, [2, 6]),
        MajoranaSum(n_sites, :nn, [1, 4]) + MajoranaSum(n_sites, :n, 7)]
        @test _circuit_matches_fr(specs, false, obs)
    end
end

@testset "sequential gates (spinful, mixed)" begin
    n_sites = 5
    specs = [
        (:hopup, [1, 2], 0.3), (:hopdn, [2, 3], -0.25), (:nupndn, [2], 0.6),
        (:hopupdn, [3, 1], 0.4), (:hopup, [4, 3], -0.5), (:pairup, [4, 5], 0.2),
        (:hole, [3], 0.35), (:hopdn, [5, 4], 0.45), (:pairdn, [1, 2], -0.3),
        (:nupndn, [4], 0.7), (:hopupdn, [2, 5], -0.15),
    ]
    for obs in [MajoranaSum(n_sites, :nup, 3),
        MajoranaSum(n_sites, :hopup, [2, 4]),
        MajoranaSum(n_sites, :nupndn, 1) + MajoranaSum(n_sites, :ndn, 5)]
        @test _circuit_matches_fr(specs, true, obs)
    end
end

@testset "sequential random circuits" begin
    Random.seed!(20240617)

    # spinless: 25 random two-site gates applied one after the other
    let n_sites = 7
        for _ in 1:5
            specs = Tuple{Symbol,Vector{Int},Float64}[]
            for _ in 1:25
                sym = rand((:hop, :nn, :pair))
                i = rand(1:n_sites)
                j = rand(1:n_sites)
                push!(specs, (sym, [i, j], randn()))
            end
            obs = MajoranaSum(n_sites, :n, rand(1:n_sites))
            @test _circuit_matches_fr(specs, false, obs; atol=1e-7)
        end
    end

    # spinful: 20 random gates (single- and two-site, incl. direction-sensitive hopupdn)
    let n_sites = 4
        for _ in 1:5
            specs = Tuple{Symbol,Vector{Int},Float64}[]
            for _ in 1:20
                pick = rand(1:5)
                if pick == 1
                    push!(specs, (:nupndn, [rand(1:n_sites)], randn()))
                elseif pick == 2
                    push!(specs, (:hole, [rand(1:n_sites)], randn()))
                elseif pick == 3
                    push!(specs, (:hopup, [rand(1:n_sites), rand(1:n_sites)], randn()))
                elseif pick == 4
                    push!(specs, (:hopdn, [rand(1:n_sites), rand(1:n_sites)], randn()))
                else
                    push!(specs, (:hopupdn, [rand(1:n_sites), rand(1:n_sites)], randn()))
                end
            end
            obs = MajoranaSum(n_sites, :nupndn, rand(1:n_sites))
            @test _circuit_matches_fr(specs, true, obs; atol=1e-7)
        end
    end
end

# ---------------------------------------------------------------------------
# Hubbard-like Trotter circuit (mirrors examples/Hubbard_1d.ipynb)
# ---------------------------------------------------------------------------

# one first-order-Trotter Hubbard layer: spin-up and spin-down hoppings on every topology edge,
# plus on-site repulsion on every site. Returns matching FermionicRotation and lookup circuits.
function _hubbard_layer(topo, n_spinful_sites, t, U, dt)
    fr = FermionicRotation[]
    look = []
    thetas = Float64[]
    for (i, j) in topo
        push!(fr, FermionicRotation(:hopup, [i, j]))
        push!(look, FermionicRotationLookup(:hopup, 2, true, [i, j], -t * dt))
        push!(thetas, -t * dt)
    end
    for (i, j) in topo
        push!(fr, FermionicRotation(:hopdn, [i, j]))
        push!(look, FermionicRotationLookup(:hopdn, 2, true, [i, j], -t * dt))
        push!(thetas, -t * dt)
    end
    for i in 1:n_spinful_sites
        push!(fr, FermionicRotation(:nupndn, i))
        push!(look, FermionicRotationLookup(:nupndn, 1, true, [i], U * dt))
        push!(thetas, U * dt)
    end
    return fr, look, thetas
end

@testset "Hubbard-like circuit" begin
    t = 1.0
    U = 4.0
    dt = 0.1
    n_layers = 3

    @testset "1D bricklayer N=$n_spinful_sites" for n_spinful_sites in (4, 5)
        topo = bricklayertopology(n_spinful_sites)
        fr, look, thetas = _hubbard_layer(topo, n_spinful_sites, t, U, dt)

        # apply the layer `n_layers` times (Trotter evolution)
        full_fr = repeat(fr, n_layers)
        full_thetas = repeat(thetas, n_layers)
        full_look = repeat(look, n_layers)

        fock = FockState(n_spinful_sites, :checkerboard, true)
        observables = [
            MajoranaSum(n_spinful_sites, :nup, 2),
            MajoranaSum(n_spinful_sites, :nupndn, 3),
            MajoranaSum(n_spinful_sites, :ndn, 1) + MajoranaSum(n_spinful_sites, :nup, n_spinful_sites),
        ]
        for obs in observables
            for min_abs_coeffs in [1.e-4, 1.e-6, 1.e-12]
                # The lookup gate applies each gate in a single shot and truncates once per gate.
                # The standard FermionicRotation truncates after every constituent Majorana rotation
                # by default (`truncate_each_mr=true` for hoppings), which is a different truncation
                # granularity and diverges under aggressive truncation. For a like-for-like
                # comparison the reference must also truncate once per gate (`truncate_each_mr=false`).
                ref = propagate(full_fr, deepcopy(obs), full_thetas; min_abs_coeff=min_abs_coeffs, truncate_each_mr=false)
                got = propagate(full_look, deepcopy(obs); min_abs_coeff=min_abs_coeffs)
                # the full propagated Majorana sums match...
                @test _msums_match(ref, got; atol=1e-8)
                # ...and so does the physical observable (overlap with the checkerboard state)
                @test isapprox(overlapwithfock(ref, fock), overlapwithfock(got, fock); atol=1e-8)
            end
        end
    end

    @testset "2D rectangle 2x2" begin
        n_spinful_sites = 4
        topo = rectangletopology(2, 2)
        fr, look, thetas = _hubbard_layer(topo, n_spinful_sites, t, U, dt)

        full_fr = repeat(fr, n_layers)
        full_thetas = repeat(thetas, n_layers)
        full_look = repeat(look, n_layers)

        obs = MajoranaSum(n_spinful_sites, :nup, 1)
        # truncate once per gate to match the lookup gate's granularity (see the 1D case above)
        ref = propagate(full_fr, deepcopy(obs), full_thetas; min_abs_coeff=1e-12, truncate_each_mr=false)
        got = propagate(full_look, deepcopy(obs); min_abs_coeff=1e-12)
        @test _msums_match(ref, got; atol=1e-8)
    end
end

# ---------------------------------------------------------------------------
# independent ground-truth cross-check against Yao (free fermions)
# ---------------------------------------------------------------------------

@testset "lookup vs Yao (free fermions)" begin
    include("yao_helpers/fermionicgates_to_yao.jl")

    n_fermions = 8
    h = 0.2
    topo = bricklayertopology(n_fermions)

    obs = MajoranaSum(n_fermions, :n, 4)

    occupied_sites = 1:2:n_fermions
    fock_state = FockState(n_fermions, occupied_sites)

    yao_psi = Yao.zero_state(n_fermions)
    state_prep = Yao.chain(n_fermions, Yao.put(site => Yao.X) for site in occupied_sites)
    Yao.apply!(yao_psi, state_prep)

    fr_circ = FermionicRotation[]
    look_circ = []
    thetas = Float64[]
    for (i, j) in topo
        push!(fr_circ, FermionicRotation(:hop, [i, j]))
        push!(look_circ, FermionicRotationLookup(:hop, 2, false, [i, j], h))
        push!(thetas, h)
    end
    yao_circ = circ_to_yao(n_fermions, fr_circ, thetas)
    yao_obs = Yao.kron(n_fermions, 4 => Yao.Z)

    for _ in 1:25
        propagate!(look_circ, obs; min_abs_coeff=-1.0)
        Yao.apply!(yao_psi, yao_circ)

        mp_res = overlapwithfock(obs, fock_state)
        yao_res = (1.0 - Yao.expect(yao_obs, yao_psi)) / 2.0
        @test abs(mp_res - yao_res) < 1e-10
    end
end

# ---------------------------------------------------------------------------
# angle-free (symbolic) lookup: build the table once, evaluate at many angles
# ---------------------------------------------------------------------------

# `table` is a SymbolicFermionicRotationLookup built once (no theta); evaluate it at `theta` and
# check the materialized gate matches both the reference FermionicRotation and the theta-baked
# FermionicRotationLookup.
function _symbolic_matches(table, symbol, is_spinful, site_inds, theta, obs; atol=1e-9)
    ref = propagate([FermionicRotation(symbol, site_inds)], deepcopy(obs), [theta]; min_abs_coeff=1e-14)

    baked = FermionicRotationLookup(symbol, length(site_inds), is_spinful, site_inds, theta)
    baked_res = propagate([baked], deepcopy(obs); min_abs_coeff=1e-14)

    got = propagate([evaluate(table, theta)], deepcopy(obs); min_abs_coeff=1e-14)

    return _msums_match(ref, got; atol=atol) && _msums_match(baked_res, got; atol=atol)
end

@testset "angle-free (symbolic) lookup" begin
    @testset "spinless" begin
        n_sites = 6
        thetas = [0.0, 0.13, 0.4, 1.0, -0.7, 2.4]
        observables = [
            MajoranaSum(n_sites, :n, 3),
            MajoranaSum(n_sites, :nn, [3, 4]),
            MajoranaSum(n_sites, :hop, [3, 5]),
            MajoranaSum(n_sites, :pair, [2, 4]),
            MajoranaSum(n_sites, :n, 2) + MajoranaSum(n_sites, :hop, [2, 5]),
        ]
        for (sym, sites) in [(:n, [3]), (:hop, [2, 3]), (:hop, [3, 5]), (:nn, [3, 4]), (:nn, [2, 5]), (:pair, [3, 4])]
            table = FermionicRotationLookup(sym, length(sites), false, sites)   # built once, no theta
            for obs in observables, th in thetas
                @test _symbolic_matches(table, sym, false, sites, th, obs)
            end
        end
    end

    @testset "spinful" begin
        n_sites = 4
        thetas = [0.0, 0.4, -0.7, 1.9]
        observables = [
            MajoranaSum(n_sites, :nup, 2),
            MajoranaSum(n_sites, :nupndn, 2),
            MajoranaSum(n_sites, :hopup, [2, 4]),
        ]
        for (sym, sites) in [(:nupndn, [2]), (:hole, [3]), (:hopup, [1, 2]), (:hopdn, [2, 3]), (:hopupdn, [3, 1])]
            table = FermionicRotationLookup(sym, length(sites), true, sites)
            for obs in observables, th in thetas
                @test _symbolic_matches(table, sym, true, sites, th, obs)
            end
        end
    end

    @testset "table is reusable across angles" begin
        # build once, evaluate at many angles; each must match the per-angle baked gate
        obs = MajoranaSum(6, :n, 3)
        table = FermionicRotationLookup(:hop, 2, false, [2, 4])
        for th in [-1.3, 0.0, 0.2, 0.75, 1.9, 3.1]
            ref = propagate([FermionicRotationLookup(:hop, 2, false, [2, 4], th)], deepcopy(obs); min_abs_coeff=1e-14)
            got = propagate([evaluate(table, th)], deepcopy(obs); min_abs_coeff=1e-14)
            @test _msums_match(ref, got)
        end
    end
end
