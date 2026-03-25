using Base.Threads

function PropagationBase.merge!(prop_cache::MultiVectorMajoranaPropagationCache; is_gaussian, kwargs...)
    if is_gaussian
        _gaussian_merge!(prop_cache; kwargs...)
    else
        _non_gaussian_merge!(prop_cache; kwargs...)
    end
end

function _gaussian_merge!(prop_cache::MultiVectorMajoranaPropagationCache; kwargs...)
    caches = prop_caches(prop_cache)
    sorted_keys = sort(collect(keys(caches)))

    # Each weight sector has an independent cache and is always merged in parallel.
    @threads for ii in eachindex(sorted_keys)
        weight_key = sorted_keys[ii]
        merge!(caches[weight_key]; kwargs...)
    end
end

function _non_gaussian_merge!(prop_cache::MultiVectorMajoranaPropagationCache; kwargs...)
    caches = prop_caches(prop_cache)
    sorted_keys = sort(collect(keys(caches)))

    # Direction +2: wave 1 (mod4 in 0,1), wave 2 (mod4 in 2,3)
    _merge_weight_direction_wave!(caches, sorted_keys, 2, (0, 1); kwargs...)
    _merge_weight_direction_wave!(caches, sorted_keys, 2, (2, 3); kwargs...)

    # Direction -2: wave 1 (mod4 in 2,3), wave 2 (mod4 in 0,1)
    _merge_weight_direction_wave!(caches, sorted_keys, -2, (2, 3); kwargs...)
    _merge_weight_direction_wave!(caches, sorted_keys, -2, (0, 1); kwargs...)
end

function _merge_weight_direction_wave!(
    caches::Dict,
    sorted_keys::Vector{Int},
    shift::Int,
    wave_residues::Tuple{Int,Int};
    kwargs...
)
    wave_origin_keys = Int[]

    # Prepare jobs serially to avoid concurrent Dict writes.
    for weight_key in sorted_keys
        residue = mod(weight_key, 4)
        if residue != wave_residues[1] && residue != wave_residues[2]
            continue
        end

        destination_weight = weight_key + shift
        origin_cache = caches[weight_key]
        if !_has_active_terms_with_weight(origin_cache, destination_weight)
            continue
        end

        get!(caches, destination_weight) do
            new_mainsum = Base.similar(mainsum(origin_cache))
            new_cache = VectorMajoranaPropagationCache(new_mainsum)
            setactivesize!(new_cache, 0)
            return new_cache
        end
        push!(wave_origin_keys, weight_key)
    end

    if isempty(wave_origin_keys)
        return
    end

    @threads for ii in eachindex(wave_origin_keys)
        weight_key = wave_origin_keys[ii]
        destination_weight = weight_key + shift
        origin_cache = caches[weight_key]
        destination_cache = caches[destination_weight]
        _move_weight_sector!(origin_cache, destination_cache, destination_weight; kwargs...)
    end
end

function _has_active_terms_with_weight(origin_cache::VectorMajoranaPropagationCache, destination_weight::Int)
    for pstr in activeterms(origin_cache)
        if get_weight(pstr) == destination_weight
            return true
        end
    end
    return false
end

function _move_weight_sector!(
    origin_cache::VectorMajoranaPropagationCache,
    destination_cache::VectorMajoranaPropagationCache,
    destination_weight::Int;
    kwargs...
)
    strings_to_be_moved(pstr) = (get_weight(pstr) == destination_weight)
    flagterms!(strings_to_be_moved, origin_cache)
    flags = activeflags(origin_cache)
    if !any(flags)
        return
    end

    destination_cache_size = activesize(destination_cache)

    flagstoindices!(origin_cache)
    n_to_be_moved = lastactiveindex(origin_cache)

    n_new = destination_cache_size + n_to_be_moved
    resize_factor = 2
    if capacity(destination_cache) < n_new
        resize!(destination_cache, n_new * resize_factor)
    end

    indices = activeindices(origin_cache)

    origin_active_terms = activeterms(origin_cache)
    origin_terms = majoranas(mainsum(origin_cache))
    origin_coeffs = coefficients(mainsum(origin_cache))

    destination_terms = majoranas(mainsum(destination_cache))
    destination_coeffs = coefficients(mainsum(destination_cache))

    for ii in eachindex(origin_active_terms)
        if flags[ii]
            term = origin_terms[ii]
            coeff = origin_coeffs[ii]

            destination_terms[destination_cache_size + indices[ii]] = term
            destination_coeffs[destination_cache_size + indices[ii]] = coeff
        end
    end
    destination_cache.active_size += n_to_be_moved

    # Merge destination cache to combine duplicates before the next wave.
    merge!(destination_cache; kwargs...)

    # Keep terms that were not moved.
    origin_cache.flags .= .!(origin_cache.flags)
    filterviaflags!(origin_cache)
end