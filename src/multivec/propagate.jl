function PropagationBase.merge!(prop_cache::MultiVectorMajoranaPropagationCache; is_gaussian, kwargs...)
    if is_gaussian
        _gassussian_merge!(prop_cache; kwargs...)
    else
        _non_gaussian_merge!(prop_cache; kwargs...)
    end
end

function _gassussian_merge!(prop_cache::MultiVectorMajoranaPropagationCache; kwargs...)
    caches = prop_caches(prop_cache)
    sorted_keys = sort(collect(keys(caches)))

    for weight_key in sorted_keys
        origin_cache = caches[weight_key]
        merge!(origin_cache; kwargs...)
    end
end

function _non_gaussian_merge!(prop_cache::MultiVectorMajoranaPropagationCache; kwargs...)
    caches = prop_caches(prop_cache)
    sorted_keys = sort(collect(keys(caches)))

    for weight_key in sorted_keys
        origin_cache = caches[weight_key]

        for destination_weight in (weight_key + 2, weight_key - 2)
            strings_to_be_moved(pstr) = (get_weight(pstr) == destination_weight)
            flagterms!(strings_to_be_moved, origin_cache)
            flags = activeflags(origin_cache)
            if !any(flags)
                continue
            end

            destination_cache = get!(caches, destination_weight) do
                new_mainsum = Base.similar(mainsum(origin_cache))
                new_cache = VectorMajoranaPropagationCache(new_mainsum)
                setactivesize!(new_cache, 0)
                return new_cache
            end

            destination_cache_size = activesize(destination_cache)

            flagstoindices!(origin_cache)
            n_to_be_moved = lastactiveindex(origin_cache)

            n_new = destination_cache_size + n_to_be_moved
            resize_factor = 2
            if capacity(destination_cache) < n_new
                #println("Resizing destination cache for weight sector $destination_weight from capacity $(capacity(destination_cache)) to $(n_new * resize_factor).")
                resize!(destination_cache, n_new * resize_factor)
            end

            indices = activeindices(origin_cache)

            origin_active_terms = activeterms(origin_cache)
            origin_terms = majoranas(mainsum(origin_cache))
            origin_coeffs = coefficients(mainsum(origin_cache))

            destination_terms = majoranas(mainsum(destination_cache))
            destination_coeffs = coefficients(mainsum(destination_cache))

            AK.foreachindex(origin_active_terms) do ii
                if flags[ii]
                    term = origin_terms[ii]
                    coeff = origin_coeffs[ii]

                    destination_terms[destination_cache_size + indices[ii]] = term
                    destination_coeffs[destination_cache_size + indices[ii]] = coeff
                end
            end
            destination_cache.active_size += n_to_be_moved

            # merge destination cache to combine duplicates
            merge!(destination_cache; kwargs...)

            # cleaning origin_cache: keep terms that were not moved
            origin_cache.flags .= .!(origin_cache.flags)
            filterviaflags!(origin_cache)
        end
    end
end