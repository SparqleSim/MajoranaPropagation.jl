using MajoranaPropagation
using PauliPropagation

using Test
using Random
Random.seed!(42)

using ProgressMeter

function random_circuit(nfermions, n_gates, n_steps, maxW_gates, msum_terms, min_abs_coeffs)
    TT = getinttype(nfermions)
    max_val = MajoranaString(nfermions, [2 * nfermions]).gammas
    circ = MajoranaRotation{TT}[]
    thetas = Float64[]

    for _ = 1:n_gates
        ms_weight = rand(2:2:maxW_gates)
        gammas = Int[]
        while length(gammas) < ms_weight
            g = rand(1:2*nfermions)
            if g ∉ gammas
                push!(gammas, g)
            end
        end
        push!(circ, MajoranaRotation(MajoranaString(nfermions, gammas)))
        push!(thetas, rand() * 2 * π)
    end

    msum = MajoranaSum(Float64, nfermions, false)
    for _ = 1:msum_terms
        ms_int = rand(TT)
        while get_weight(ms_int) % 2 != 0 || ms_int > max_val
            ms_int = rand(TT)
        end
        PauliPropagation.PropagationBase.add!(msum, ms_int, randn())
    end

    for min_abs_coeff in min_abs_coeffs
        obs = deepcopy(msum)
        obs_vec = VectorMajoranaSum(deepcopy(msum))
        truncate!(obs; min_abs_coeff)
        truncate!(obs_vec; min_abs_coeff)
        @showprogress for _ = 1:n_steps
            @time propagate!(circ, obs, thetas; min_abs_coeff)
            @time propagate!(circ, obs_vec, thetas; min_abs_coeff)
            @test length(obs) == length(obs_vec)
            @test obs == obs_vec
        end
    end
end



@testset "Vector propagation - Hubbard" begin
    N_sites = 20
    topo = bricklayertopology(N_sites)

    U = 1.
    t = 1.
    dτ = 0.1
    n_steps = 12

    fock_state = FockState(N_sites, :checkerboard, true) #create a checkerboard state with spinful fermions

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

    for min_abs_coeff in min_abs_coeffs
        obs = MajoranaSum(N_sites, :nupndn, 8)
        MajoranaPropagation.pop_id!(obs)
        obs_vec = VectorMajoranaSum(deepcopy(obs))
        for _ = 1:n_steps
            propagate!(circ, obs, thetas; min_abs_coeff)
            propagate!(circ, obs_vec, thetas; min_abs_coeff)
            @test length(obs) == length(obs_vec)
            @test obs == obs_vec
            @test abs(overlapwithfock(obs, fock_state) - overlapwithfock(obs_vec, fock_state)) < 1.e-14
        end
    end
end