using MajoranaPropagation
using PauliPropagation
using Base.Threads

using Test
using Random
Random.seed!(42)

@testset "MultiDict propagation - Hubbard" begin
    N_sites = 31
    topo = bricklayertopology(N_sites)

    U = 2.
    t = 1.
    dτ = 0.12
    n_steps = 10

    fock_state = FockState(N_sites, 2:2:N_sites, 1:2:N_sites)

    circ = FermionicGate[]
    thetas = Float64[]

    for pair in topo
        push!(circ, FermionicGate(:hopup, pair))
        push!(thetas, -t * dτ)
    end
    for pair in topo
        push!(circ, FermionicGate(:hopdn, pair))
        push!(thetas, -t * dτ)
    end
    for i = 1:N_sites
        push!(circ, FermionicGate(:nupndn, i))
        push!(thetas, U * dτ)
    end

    min_abs_coeffs = [1.e-2, 1.e-3, 1.e-4, 1.e-5]

    n_levels = nthreads()
    level_mapper = ms -> mod(mod(ms * 11, n_levels) + countweight(ms), n_levels) + 1


    for min_abs_coeff in min_abs_coeffs
        obs = MajoranaSum(N_sites, :ndn, 21)
        MajoranaPropagation.pop_id!(obs)
        obs_multidict = MajoranaSumMulti(deepcopy(obs))
        for _ = 1:n_steps
            propagate!(circ, obs, thetas; min_abs_coeff)
            propagate!(circ, obs_multidict, thetas; min_abs_coeff)
            @test length(obs) == length(obs_multidict)
            @test obs == obs_multidict
            @test abs(overlapwithfock(obs, fock_state) - overlapwithfock(obs_multidict, fock_state)) < 1.e-14
        end
    end
end