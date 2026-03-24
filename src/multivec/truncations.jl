function PropagationBase.truncate!(prop_cache::MultiVectorMajoranaPropagationCache; kwargs...)
    for (weight_sector, vpropcache) in prop_caches(prop_cache)
        #=if weight_sector == 2
            println("---- in truncations")
            @show length(vpropcache)
        end=#
        truncate!(vpropcache; kwargs...)
        #=if weight_sector == 2
            @show length(vpropcache)
        end=#
    end
    return prop_cache
end