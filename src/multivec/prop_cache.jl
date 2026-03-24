
# =========================
# Propagation cache
# =========================

mutable struct MultiVectorMajoranaPropagationCache{VMS<:VectorMajoranaSum,VB,VI} <: AbstractMajoranaPropagationCache
    prop_caches::Dict{Int,VectorMajoranaPropagationCache{VMS,VB,VI}}
end

prop_caches(prop_cache::MultiVectorMajoranaPropagationCache) = prop_cache.prop_caches

function PropagationBase.mainsum(prop_cache::MultiVectorMajoranaPropagationCache{VectorMajoranaSum{Vector{TT},Vector{CT}},VB,VI}) where {VB,VI, TT<:Integer,CT}
    mouts = Dict{Int,VectorMajoranaSum{Vector{TT},Vector{CT}}}()
    for (weight_sector, vpropcache) in prop_caches(prop_cache)
        mouts[weight_sector] = mainsum(vpropcache)
    end
    sample_sum = first(values(mouts))
    return MultiVectorMajoranaSum(PropagationBase.nsites(sample_sum), is_spinful(sample_sum), mouts)
end
#PropagationBase.auxsum(prop_cache::MultiVectorMajoranaPropagationCache) = prop_cache.aux_msum_dict
#nfermions(prop_cache::MultiVectorMajoranaPropagationCache) = nfermions(mainsum(prop_cache))

#=function PropagationBase.setmainsum!(prop_cache::MultiVectorMajoranaPropagationCache, msum::MultiVectorMajoranaSum)
    prop_cache.main_msum = msum
    return prop_cache
end

function PropagationBase.setauxsum!(prop_cache::MultiVectorMajoranaPropagationCache, aux_msum::MultiVectorMajoranaSum)
    prop_cache.aux_msum = aux_msum
    return prop_cache
end=#

"""
    MultiVectorMajoranaPropagationCache(msum::MultiVectorMajoranaSum)

Create a propagation cache for a `MultiVectorMajoranaSum`. For each
weight sector in the main sum we allocate an auxiliary multi-vector sum
with support in the neighbouring weight sectors.
"""
function MultiVectorMajoranaPropagationCache(multimsum::MultiVectorMajoranaSum{TT,CT}) where {TT<:Integer,CT}
    prop_caches = Dict{Int,VectorMajoranaPropagationCache{VectorMajoranaSum{Vector{TT},Vector{CT}},Vector{Bool},Vector{Int}}}()

    for (w, vms) in multimsum.MultiMajoranas
        prop_caches[w] = VectorMajoranaPropagationCache(vms)
    end
    return MultiVectorMajoranaPropagationCache(prop_caches)
end

PropagationBase.PropagationCache(multimsum::MultiVectorMajoranaSum) = MultiVectorMajoranaPropagationCache(multimsum)

function nsites(prop_cache::MultiVectorMajoranaPropagationCache)
    #pd = mainsum(prop_cache)[first(keys(mainsum(prop_cache)))]
    #@show typeof(pd)
    #@show nsites(pd)
    #@show is_spinful(pd)
    #@show nfermions(pd)
    return PropagationBase.nsites(prop_caches(prop_cache)[first(keys(prop_caches(prop_cache)))])
end

function Base.length(prop_cache::MultiVectorMajoranaPropagationCache)
    total_length = 0
    for (weight_sector, vpropcache) in prop_caches(prop_cache)
        total_length += activesize(vpropcache)
    end
    return total_length
end

function show_stats(prop_cache::MultiVectorMajoranaPropagationCache)
    sorted_keys = sort(collect(keys(prop_caches(prop_cache))))
    println("MultiVectorMajoranaPropagationCache stats:")
    tot_strings = 0 
    for weight_sector in sorted_keys
        vpropcache = prop_caches(prop_cache)[weight_sector]
        println("  Weight sector $weight_sector: $(activesize(vpropcache)) strings")
        tot_strings += activesize(vpropcache)
    end
    println("Total strings: $tot_strings")
end