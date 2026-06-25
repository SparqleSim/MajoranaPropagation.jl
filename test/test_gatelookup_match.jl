# GateLookup vs standard propagation.
#
# For a range of fermionic gates, both lookup variants must reproduce the standard
# `FermionicRotation` propagation:
#   - the theta-baked `FermionicRotationLookup(symbol, …, theta)` (a StaticGate), and
#   - the angle-free `FermionicRotationLookup(symbol, …)` built once and `evaluate`d at `theta`.
#
# `truncate_each_mr=false` is set everywhere. A `FermionicRotationLookup` applies the whole gate in
# one shot, so to compare like-for-like the reference `FermionicRotation` must NOT truncate between
# its constituent Majorana rotations. This matters for gates whose default is `truncate_each_mr=true`
# (everything except :hop/:hopup/:hopdn, see `flag_non_number_preserving`). `min_abs_coeff=-1.0`
# disables coefficient truncation so the comparison isolates the gate action itself.

using MajoranaPropagation
using MajoranaPropagation.GateLookup
using Test

# coefficient dict keyed by BigInt so sums with different integer widths compare cleanly
function _gl_coeffdict(msum)
    d = Dict{BigInt,ComplexF64}()
    for (k, v) in msum
        d[BigInt(k)] = ComplexF64(v)
    end
    return d
end

function _gl_match(a, b; atol=1e-9)
    da, db = _gl_coeffdict(a), _gl_coeffdict(b)
    for k in union(keys(da), keys(db))
        if abs(get(da, k, 0.0im) - get(db, k, 0.0im)) > atol
            return false
        end
    end
    return true
end

# `sym_table` is the angle-free lookup, built once per gate and reused across angles. Propagate
# `obs` through the standard FermionicRotation, the theta-baked lookup, and the evaluated symbolic
# lookup, and assert all three agree.
function _check_lookup(sym_table, symbol, is_spinful, sites, theta, obs; atol=1e-9)
    ref = propagate([FermionicRotation(symbol, sites)], deepcopy(obs), [theta];
        truncate_each_mr=false, min_abs_coeff=-1.0)

    baked = FermionicRotationLookup(symbol, length(sites), is_spinful, sites, theta)
    baked_res = propagate([baked], deepcopy(obs); truncate_each_mr=false, min_abs_coeff=-1.0)

    sym_res = propagate([evaluate(sym_table, theta)], deepcopy(obs);
        truncate_each_mr=false, min_abs_coeff=-1.0)

    @test _gl_match(ref, baked_res; atol=atol)   # theta-baked lookup == standard
    @test _gl_match(ref, sym_res; atol=atol)     # angle-free lookup == standard
end

@testset "GateLookup matches standard propagation (truncate_each_mr=false)" begin
    thetas = [0.0, 0.3, -0.8, 1.7]

    @testset "spinless gates" begin
        n = 6
        observables = [
            MajoranaSum(n, :n, 3),
            MajoranaSum(n, :hop, [2, 5]),                            # spans the gate boundary
            MajoranaSum(n, :nn, [3, 4]) + MajoranaSum(n, :n, 1),
        ]
        gates = [
            (:hop, [3, 4]),
            (:hop, [2, 5]),     # non-contiguous
            (:nn, [3, 4]),
            (:nn, [2, 5]),
            (:pair, [3, 4]),
        ]
        for (sym, sites) in gates
            table = FermionicRotationLookup(sym, length(sites), false, sites)  # built once
            for obs in observables, th in thetas
                _check_lookup(table, sym, false, sites, th, obs)
            end
        end
    end

    @testset "spinful gates" begin
        n = 4
        observables = [
            MajoranaSum(n, :nup, 2),
            MajoranaSum(n, :nupndn, 2),
            MajoranaSum(n, :hopup, [1, 3]) + MajoranaSum(n, :ndn, 4),
        ]
        gates = [
            (:hopup, [1, 2]),
            (:hopdn, [2, 3]),
            (:hopupdn, [3, 1]),   # direction-sensitive, reversed site order
            (:nupndn, [2]),
            (:hole, [3]),
            (:pairup, [1, 2]),
            (:pairdn, [2, 3]),
        ]
        for (sym, sites) in gates
            table = FermionicRotationLookup(sym, length(sites), true, sites)  # built once
            for obs in observables, th in thetas
                _check_lookup(table, sym, true, sites, th, obs)
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Canonical (placement-free) lookup: build once per gate type, place on any sites
# ---------------------------------------------------------------------------

# `table` is a CanonicalFermionicRotationLookup built once for a gate type. Place it on `site_inds`
# at `theta` and assert it matches both the standard FermionicRotation and the site_inds-baked lookup.
function _check_canonical(table, symbol, is_spinful, site_inds, theta, obs; atol=1e-9)
    ref = propagate([FermionicRotation(symbol, site_inds)], deepcopy(obs), [theta];
        truncate_each_mr=false, min_abs_coeff=-1.0)

    baked = FermionicRotationLookup(symbol, length(site_inds), is_spinful, site_inds, theta)
    baked_res = propagate([baked], deepcopy(obs); truncate_each_mr=false, min_abs_coeff=-1.0)

    canon = propagate([evaluate(table, theta, site_inds)], deepcopy(obs);
        truncate_each_mr=false, min_abs_coeff=-1.0)

    @test _gl_match(ref, canon; atol=atol)        # canonical placement == standard
    @test _gl_match(baked_res, canon; atol=atol)  # canonical placement == site_inds-baked lookup
end

@testset "Canonical (placement-free) lookup" begin
    thetas = [0.0, 0.35, -0.9, 1.6]

    @testset "spinless" begin
        n = 6
        observables = [
            MajoranaSum(n, :n, 3),
            MajoranaSum(n, :hop, [2, 5]),
            MajoranaSum(n, :nn, [3, 4]) + MajoranaSum(n, :n, 1),
        ]
        for (sym, placements) in [
            (:hop, [[1, 2], [3, 5], [5, 3], [2, 6]]),   # incl. reversed order
            (:nn, [[2, 3], [4, 1]]),
            (:pair, [[1, 4], [4, 1]]),
        ]
            table = FermionicRotationLookup(sym, 2, false)   # canonical, built once per gate type
            for site_inds in placements, obs in observables, th in thetas
                _check_canonical(table, sym, false, site_inds, th, obs)
            end
        end
    end

    @testset "spinful" begin
        n = 4
        observables = [
            MajoranaSum(n, :nup, 2),
            MajoranaSum(n, :nupndn, 3),
            MajoranaSum(n, :hopup, [1, 3]),
        ]
        for (sym, placements) in [
            (:hopup, [[1, 2], [2, 4], [4, 2]]),     # symmetric: any order
            (:hopdn, [[1, 3], [3, 1]]),
            (:hopupdn, [[1, 2], [1, 3], [2, 4]]),   # direction-sensitive: ascending placement
            (:pairup, [[1, 2], [2, 4]]),
        ]
            table = FermionicRotationLookup(sym, 2, true)
            for site_inds in placements, obs in observables, th in thetas
                _check_canonical(table, sym, true, site_inds, th, obs)
            end
        end
        for sym in (:nupndn, :hole)   # single-site gates
            table = FermionicRotationLookup(sym, 1, true)
            for site in [1, 2, 4], obs in observables, th in thetas
                _check_canonical(table, sym, true, [site], th, obs)
            end
        end
    end

    @testset "invalid site_inds error" begin
        table = FermionicRotationLookup(:hop, 2, false)
        @test_throws ArgumentError evaluate(table, 0.3, [2, 2])      # repeated indices
        @test_throws ArgumentError evaluate(table, 0.3, [2])         # too few sites
        @test_throws ArgumentError evaluate(table, 0.3, [1, 2, 3])   # too many sites
    end
end
