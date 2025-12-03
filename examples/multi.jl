using Revise
using MajoranaPropagation
using PauliPropagation

using Plots 
using LaTeXStrings
using ProgressMeter
using TimerOutputs
using Base.Threads

@show nthreads()


let 
    nx = 7 
    ny = 7
    nspinful = nx * ny 
    topo = rectangletopology(nx, ny)

    U = 1.
    t = 2.
    dt = 0.06 

    circ_single = []
    thetas_single = []

    #up hoppings 
    for (i, j) in topo
        push!(circ_single, FermionicGate(:hopup, [i, j]))
        push!(thetas_single, -t * dt)
    end

    #down hoppings 
    for (i, j) in topo
        push!(circ_single, FermionicGate(:hopdn, [i, j]))
        push!(thetas_single, -t * dt)
    end

    #on-site repulsion 
    for i = 1:nspinful
        push!(circ_single, FermionicGate(:nupndn, i))
        push!(thetas_single, U * dt)
    end

    msum = MajoranaSum(nspinful, :nupndn, 3)
    id_val = pop_id!(msum)
    multi_msum = MajoranaSumMulti(msum)

    @show multi_msum.MultiMajoranas

    min_abs_coeff = 1.e-6

    n_reps = 5

    times_multi = zeros(n_reps)
    times_normal = zeros(n_reps)
    lengths_multi = zeros(n_reps)
    lengths_normal = zeros(n_reps)

    for k = 1:n_reps
        println("---$(k)---")
        # multisum 
        times_multi[k] = @elapsed  propagate!(circ_single, multi_msum, thetas_single; min_abs_coeff=min_abs_coeff)
        #println(gfhj)
        #normal mode 
        times_normal[k] = @elapsed  propagate!(circ_single, msum, thetas_single; min_abs_coeff=min_abs_coeff)
        @show times_multi[k]
        @show times_normal[k]
        show_stats(multi_msum)
        #@show length(multi_msum)
        @show length(msum)
        lengths_multi[k] = length(multi_msum)
        lengths_normal[k] = length(msum)
    end

    #profile 
    #@profview propagate!(circ_single, multi_msum, thetas_single; min_abs_coeff=min_abs_coeff)

end