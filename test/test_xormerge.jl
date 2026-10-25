using MajoranaPropagation
using PauliPropagation
using Test
using Random

# End-to-end checks of the XOR merge that Majorana rotations go through (PropagationBase's
# `xorbranchmergeandtruncate!`): the dict and vector backends must agree. The merge only spawns tasks
# when Julia runs with more than one thread; run the suite with JULIA_NUM_THREADS>1 to exercise it.

const xorm_rng = MersenneTwister(42)

# O(N) agreement check (the AbstractTermSum `==` re-merges copies and looks terms up
# one by one, which is far too slow at the term counts these blocks reach)
function xorm_sums_agree(msum, vmsum; atol=1e-10)
    d_dict = Dict(k => v for (k, v) in msum)
    d_vec = Dict(k => v for (k, v) in vmsum)
    return length(d_dict) == length(d_vec) &&
           all(isapprox(v, get(d_vec, k, Inf); atol) for (k, v) in d_dict)
end

@testset "XOR merge end-to-end (dict vs vector)" begin
    # interacting 2D problem large enough that per-rotation tails exceed
    # 2 * _MIN_ELEMS_PER_TASK, so threaded runs merge with multiple tasks
    nx, ny = 5, 4
    topo = rectangletopology(nx, ny)
    nf = nx * ny
    circ = FermionicRotation[]
    thetas = Float64[]
    for pair in topo
        push!(circ, FermionicRotation(:hop, pair))
        push!(thetas, 0.1)
        push!(circ, FermionicRotation(:nn, pair))
        push!(thetas, 0.1)
    end

    obs = MajoranaSum(nf, :n, (nf + 1) ÷ 2)
    obs_vec = VectorMajoranaSum(deepcopy(obs))
    for _ in 1:4
        propagate!(circ, obs, thetas; min_abs_coeff=1e-6)
        propagate!(circ, obs_vec, thetas; min_abs_coeff=1e-6)
        @test length(obs) == length(obs_vec)
        @test xorm_sums_agree(obs, obs_vec)
    end
    @test length(obs_vec) > 100_000

    # bare MajoranaRotation circuit: exercises the applymergetruncate! overload that
    # routes bare rotations through the gate-aware merge
    nf = 10
    TT = getinttype(nf)
    mr_circ = MajoranaRotation{TT}[]
    mr_thetas = Float64[]
    for _ in 1:40
        w = rand(xorm_rng, (2, 4))
        gammas = Int[]
        while length(gammas) < w
            gg = rand(xorm_rng, 1:2*nf)
            gg in gammas || push!(gammas, gg)
        end
        push!(mr_circ, MajoranaRotation(MajoranaString(nf, gammas)))
        push!(mr_thetas, randn(xorm_rng))
    end
    obs = MajoranaSum(nf, :n, 5)
    obs_vec = VectorMajoranaSum(deepcopy(obs))
    for _ in 1:2
        propagate!(mr_circ, obs, mr_thetas; min_abs_coeff=1e-6)
        propagate!(mr_circ, obs_vec, mr_thetas; min_abs_coeff=1e-6)
        @test length(obs) == length(obs_vec)
        @test xorm_sums_agree(obs, obs_vec)
    end
end
