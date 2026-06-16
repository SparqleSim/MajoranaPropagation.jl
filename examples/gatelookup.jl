using Revise
using MajoranaPropagation
using MajoranaPropagation.GateLookup

# `FermionicRotationLookup` is a `StaticGate` alternative to `FermionicRotation`: the rotation
# angle is baked into a precomputed transfer map (lookup table) at construction time, and the
# gate's conjugate action is then read off the table instead of decomposing into Majorana
# rotations on every application.
#
# Constructor: FermionicRotationLookup(symbol, sites_acted_on, is_spinful, site_inds, theta)

let
    # -------- spinless example --------
    n_sites = 5
    theta = 0.4

    obs = MajoranaSum(n_sites, :hop, [1,2])

    # reference: the usual parametrized FermionicRotation
    @time ref = propagate([FermionicRotation(:nn, [2, 3])], deepcopy(obs), [theta])

    # lookup-table gate (StaticGate => no thetas passed to `propagate`)
    look_gate = FermionicRotationLookup(:nn, 2, false, [2, 3], theta)
    @time got = propagate([look_gate], deepcopy(obs))

    #@show ref
    #@show got
    @assert ref == got

    # -------- spinful example --------
    n_sites = 4
    obs = MajoranaSum(n_sites, :nup, 1)

    @time ref = propagate([FermionicRotation(:hopup, [1, 2])], deepcopy(obs), [theta])
    @time got = propagate([FermionicRotationLookup(:hopup, 2, true, [1, 2], theta)], deepcopy(obs))

    #@show ref
    #@show got
    @assert ref == got 
end
