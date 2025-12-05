using Revise
using MajoranaPropagation
using PauliPropagation

#using Plots 
#using LaTeXStrings
#using ProgressMeter
#using TimerOutputs
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
    nx = 6 
    ny = 6
    nspinful = nx * ny 
    topo = rectangletopology(nx, ny)
    #nspinful = 40
    #topo = bricklayertopology(nspinful)

    U = 4.
    t = 1.
    dt = 0.06 

    circ_single = []
    thetas_single = []

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

    msum = MajoranaSum(nspinful, :nupndn, 3) #* MajoranaSum(nspinful, :nupndn, 5)
    @show msum 
    id_val = pop_id!(msum)
    multi_msum = MajoranaSumMulti(msum)

    @show multi_msum.MultiMajoranas

    min_abs_coeff = 5.e-7
    max_single_filter = create_max_single_filter(2 * nspinful)
    max_singles = 8

    if max_singles < Inf 
        custom_trunc = let max_single_filter=max_single_filter, max_singles = max_singles
            (mstr, coeff) -> (compute_max_single(mstr, 0, max_single_filter) > max_singles)
        end
    else
        custom_trunc = nothing
    end

    n_reps = 10

    times_multi = zeros(n_reps)
    times_normal = zeros(n_reps)
    lengths_multi = zeros(n_reps)
    lengths_normal = zeros(n_reps)

    for k = 1:n_reps
        println("---$(k)---")
        # multisum 
        times_multi[k] = @elapsed  propagate!(circ_single, multi_msum, thetas_single; min_abs_coeff=min_abs_coeff, customtruncfunc=custom_trunc)
        println("time multi: $(print_time(times_multi[k]))")
        #println(gfhj)
        #normal mode 
        times_normal[k] = @elapsed  propagate!(circ_single, msum, thetas_single; min_abs_coeff=min_abs_coeff, customtruncfunc=custom_trunc)
        println("time normal: $(print_time(times_normal[k]))")
        show_stats(multi_msum)
        #@show length(multi_msum)
        @show length(msum)
        lengths_multi[k] = length(multi_msum)
        lengths_normal[k] = length(msum)
    end

    #profile 
    #@profview propagate!(circ_single, multi_msum, thetas_single; min_abs_coeff=min_abs_coeff, customtruncfunc=custom_trunc)

end