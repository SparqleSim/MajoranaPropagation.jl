using Pkg
Pkg.activate("./")

using MajoranaPropagation
using PauliPropagation
using Base.Threads

@show nthreads()

function print_time(seconds)
    hours = div(seconds, 3600)
    minutes = div(seconds % 3600, 60)
    seconds = seconds % 60
    if hours > 0
        return "$(round(Int, hours))h $(round(Int, minutes))m $(round(Int, seconds))s"
    elseif minutes > 0
        return "$(round(Int, minutes))m $(round(Int, seconds))s"
    else
        return "$(seconds)s"
    end
end

let
    nx = 7
    ny = 7
    nspinful = nx * ny
    topo = rectangletopology(nx, ny)

    U = 8.0
    t = 1.0
    dt = 0.06

    circ_single = []
    thetas_single = []

    # up hoppings
    for (i, j) in topo
        push!(circ_single, FermionicGate(:hopup, [i, j]))
        push!(thetas_single, -t * dt / 2.0)
    end

    # down hoppings
    for (i, j) in topo
        push!(circ_single, FermionicGate(:hopdn, [i, j]))
        push!(thetas_single, -t * dt / 2.0)
    end

    # on-site repulsion
    for i = 1:nspinful
        push!(circ_single, FermionicGate(:nupndn, i))
        push!(thetas_single, U * dt)
    end

    # down hoppings
    for (i, j) in reverse(topo)
        push!(circ_single, FermionicGate(:hopdn, [i, j]))
        push!(thetas_single, -t * dt / 2.0)
    end

    # up hoppings
    for (i, j) in reverse(topo)
        push!(circ_single, FermionicGate(:hopup, [i, j]))
        push!(thetas_single, -t * dt / 2.0)
    end

    # initial observable
    msum = MajoranaSum(nspinful, :nupndn, 3)
    id_val = MajoranaPropagation.pop_id!(msum)

    multi_msum = MajoranaSumMulti(msum)
    vec_msum = VectorMajoranaSum(msum)
    multivec_msum = MultiVectorMajoranaSum(msum)

    min_abs_coeff = 5e-7
    max_singles = 8

    n_reps = 5

    times_multi = zeros(n_reps)
    times_vec = zeros(n_reps)
    times_multivec = zeros(n_reps)

    lengths_multi = zeros(n_reps)
    lengths_vec = zeros(n_reps)
    lengths_multivec = zeros(n_reps)

    for k = 1:n_reps
        println("---$(k)---")

        # multi dict
        #times_multi[k] = @elapsed propagate!(circ_single, multi_msum, thetas_single; min_abs_coeff=min_abs_coeff, max_unpaired=max_singles)
        #println("time multi: $(print_time(times_multi[k]))")

        # vector
        times_vec[k] = @elapsed vec_msum = propagate!(circ_single, vec_msum, thetas_single; min_abs_coeff=min_abs_coeff, max_unpaired=max_singles)
        println("time vec: $(print_time(times_vec[k]))")

        # multi-vector
        times_multivec[k] = @elapsed propagate!(circ_single, multivec_msum, thetas_single; min_abs_coeff=min_abs_coeff, max_unpaired=max_singles)
        println("time multivec: $(print_time(times_multivec[k]))")

        println("stats multi:")
        show_stats(multi_msum)
        println("stats multivec:")
        show_stats(multivec_msum)

        @show length(multi_msum)
        @show length(vec_msum)
        @show length(multivec_msum)

        #@assert length(multivec_msum) == length(multi_msum)
        @assert length(multivec_msum) == length(vec_msum)

        lengths_multi[k] = length(multi_msum)
        lengths_vec[k] = length(vec_msum)
        lengths_multivec[k] = length(multivec_msum)
    end
end

