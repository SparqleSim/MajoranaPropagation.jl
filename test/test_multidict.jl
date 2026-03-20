using Test
using MajoranaPropagation
using PauliPropagation

using TimerOutputs

@testset "MajoranaSumMulti correctness" begin
    nx = 4
    ny = 4
    nspinful = nx * ny
    topo = rectangletopology(nx, ny)

    U = 8.
    t = 1.
    dt = 0.06

    circ_single = FermionicGate[]
    thetas_single = Float64[]

    #up hoppings 
    for (i, j) in topo
        push!(circ_single, FermionicGate(:hopup, [i, j]))
        push!(thetas_single, -t * dt / 2.)
    end

    #down hoppings 
    for (i, j) in topo
        push!(circ_single, FermionicGate(:hopdn, [i, j]))
        push!(thetas_single, -t * dt / 2.)
    end

    #on-site repulsion 
    for i = 1:nspinful
        push!(circ_single, FermionicGate(:nupndn, i))
        push!(thetas_single, U * dt)
    end

    #down hoppings 
    for (i, j) in reverse(topo)
        push!(circ_single, FermionicGate(:hopdn, [i, j]))
        push!(thetas_single, -t * dt / 2.)
    end

    #up hoppings 
    for (i, j) in reverse(topo)
        push!(circ_single, FermionicGate(:hopup, [i, j]))
        push!(thetas_single, -t * dt / 2.)
    end

    min_abs_coeff = 1.e-7
    unpaired_mask = create_unpaired_mask(2 * nspinful)
    max_unpaired = 12
    n_levels = div(max_unpaired, 2) + 1
    level_mapper(mstr::TT) where {TT<:Integer} = min(div(compute_unpaired(mstr, unpaired_mask), 2) + 1, n_levels)

    msum = MajoranaSum(nspinful, :nupndn, 3) #* MajoranaSum(nspinful, :nupndn, 5)
    id_val = MajoranaPropagation.pop_id!(msum)
    multi_msum = MajoranaSumMulti(msum, level_mapper, n_levels)
    multi_sum_prop_cache = MajoranaMultiPropagationCache(multi_msum, n_levels)
    to = TimerOutput()

    n_reps = 4

    for k = 1:n_reps
        propagate!(circ_single, multi_sum_prop_cache, thetas_single; min_abs_coeff, max_unpaired, unpaired_mask, level_mapper, n_levels, to)

        propagate!(circ_single, msum, thetas_single; min_abs_coeff, max_unpaired, unpaired_mask, to)
        @test length(mainsum(multi_sum_prop_cache)) == length(msum)
    end
end