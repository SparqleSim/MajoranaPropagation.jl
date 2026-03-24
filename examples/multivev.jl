using Pkg
Pkg.activate("./")

using Revise

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
nx = 3
ny = 3
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
msum = MajoranaSum(nspinful, :nupndn, 3) * MajoranaSum(nspinful, :nupndn, 4) #* MajoranaSum(nspinful, :nupndn, 1) * MajoranaSum(nspinful, :hopup, [5, 6])
#msum = MajoranaSum(nspinful, :nup, 3)
id_val = MajoranaPropagation.pop_id!(msum)

multi_msum = MajoranaSumMulti(msum)
vec_msum = VectorMajoranaSum(msum)
multivec_msum = MultiVectorMajoranaSum(msum)
multivec_msum = MultiVectorMajoranaPropagationCache(multivec_msum)

#circ_single = [FermionicGate(:hopup, [1, 3])]
#thetas_single = [0.3]
#push!(circ_single, FermionicGate(:hopup, [3, 5]))
#push!(thetas_single, -t * dt / 2.0)


min_abs_coeff = 1.e-6
max_singles = 8

n_reps = 6

times_multi = zeros(n_reps)
times_vec = zeros(n_reps)
times_multivec = zeros(n_reps)

lengths_multi = zeros(n_reps)
lengths_vec = zeros(n_reps)
lengths_multivec = zeros(n_reps)

for k = 1:n_reps
    println("---$(k)---")

    # normal 
    propagate!(circ_single, msum, thetas_single; min_abs_coeff=min_abs_coeff, max_unpaired=max_singles)

    # vector
    times_vec[k] = @elapsed propagate!(circ_single, vec_msum, thetas_single; min_abs_coeff=min_abs_coeff, max_unpaired=max_singles)
    #println("time vec: $(print_time(times_vec[k]))")

    # multi-vector
    times_multivec[k] = @elapsed propagate!(circ_single, multivec_msum, thetas_single; min_abs_coeff=min_abs_coeff, max_unpaired=max_singles)
    #println("time multivec: $(print_time(times_multivec[k]))")

    #println("stats multi:")
    #show_stats(multi_msum)
    #println("stats multivec:")
    println("-----")
    show_stats(MajoranaSumMulti(msum))
    println("-----")
    show_stats(multivec_msum)


    @show length(msum)
    @show length(multivec_msum)
    @show length(vec_msum)

    @assert length(vec_msum) == length(multivec_msum)
end
end

