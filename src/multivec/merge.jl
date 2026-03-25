# No-AK merge path for vector caches used by MultiVectorMajoranaPropagationCache.
# This mirrors the vectorbackend merge stages but uses custom ThreadPools pools.

function _merge_vector_cache_noak!(prop_cache::VectorMajoranaPropagationCache; pool)
    n_active = activesize(prop_cache)
    if n_active <= 1
        return prop_cache
    end

    _sortbyterm_pool!(prop_cache; pool)
    _deduplicate_pool!(prop_cache; pool)

    return prop_cache
end

function _sortbyterm_pool!(prop_cache::VectorMajoranaPropagationCache; pool)
    idx_view = activeindices(prop_cache)
    term_view = activeterms(prop_cache)

    tforeach(pool, eachindex(idx_view)) do ii
        idx_view[ii] = ii
    end

    sortperm!(idx_view, term_view)
    _permuteviaindices_pool!(prop_cache; pool)

    return prop_cache
end

function _permuteviaindices_pool!(prop_cache::VectorMajoranaPropagationCache; pool)
    indices_view = activeindices(prop_cache)
    term_view = activeterms(prop_cache)
    coeffs_view = activecoeffs(prop_cache)
    aux_terms_view = activeauxterms(prop_cache)
    aux_coeffs_view = activeauxcoeffs(prop_cache)

    tforeach(pool, eachindex(indices_view)) do ii
        sorted_idx = indices_view[ii]
        aux_terms_view[ii] = term_view[sorted_idx]
        aux_coeffs_view[ii] = coeffs_view[sorted_idx]
    end

    swapsums!(prop_cache)

    return prop_cache
end

function _deduplicate_pool!(prop_cache::VectorMajoranaPropagationCache; pool)
    _flaggroupbegin_pool!(prop_cache; pool)
    group_count = _flagstoindices_pool!(prop_cache; pool)
    _mergegroups_pool!(prop_cache; pool, group_count)
    return prop_cache
end

function _flaggroupbegin_pool!(prop_cache::VectorMajoranaPropagationCache; pool)
    term_view = activeterms(prop_cache)
    flags_view = activeflags(prop_cache)

    tforeach(pool, eachindex(term_view)) do ii
        if ii == 1
            flags_view[ii] = true
        else
            flags_view[ii] = term_view[ii] != term_view[ii - 1]
        end
    end

    return prop_cache
end

function _flagstoindices_pool!(prop_cache::VectorMajoranaPropagationCache; pool)
    dst_indices = activeindices(prop_cache)
    flags = activeflags(prop_cache)
    group_starts = prop_cache.group_starts

    n = length(flags)
    n == 0 && return 0

    nblocks = min(n, 8 * Threads.nthreads())
    block_size = cld(n, nblocks)

    block_sums = prop_cache.block_sums
    block_offsets = prop_cache.block_offsets

    # Pass 1: per-block counts in parallel.
    tforeach(pool, 1:nblocks) do b
        lo = (b - 1) * block_size + 1
        hi = min(b * block_size, n)
        local_count = 0
        for ii in lo:hi
            local_count += flags[ii] ? 1 : 0
        end
        block_sums[b] = local_count
    end

    # Pass 2: sequential scan over block sums.
    running = 0
    for b in 1:nblocks
        block_offsets[b] = running
        running += block_sums[b]
    end

    # Pass 3: per-block local scans in parallel.
    # Also fill `group_starts[group_id] = start_index`.
    tforeach(pool, 1:nblocks) do b
        lo = (b - 1) * block_size + 1
        hi = min(b * block_size, n)
        local_count = block_offsets[b]
        for ii in lo:hi
            if flags[ii]
                local_count += 1
                dst_indices[ii] = local_count
                group_starts[local_count] = ii
            else
                dst_indices[ii] = local_count
            end
        end
    end

    return running
end

function _mergegroups_pool!(prop_cache::VectorMajoranaPropagationCache; pool, group_count::Int)
    term_view = activeterms(prop_cache)
    coeffs = activecoeffs(prop_cache)
    aux_terms = activeauxterms(prop_cache)
    aux_coeffs = activeauxcoeffs(prop_cache)
    active_size = activesize(prop_cache)
    group_starts = prop_cache.group_starts

    group_count <= 0 && return prop_cache

    tforeach(pool, 1:group_count) do kk
        ii = group_starts[kk]
        end_idx = kk < group_count ? group_starts[kk + 1] - 1 : active_size

        CT = typeof(coeffs[ii])
        merged_coeff = zero(CT)
        for jj in ii:end_idx
            merged_coeff = mergefunc(merged_coeff, coeffs[jj])
        end

        # Group id `kk` corresponds to the position in `aux_*` after dedup.
        aux_terms[kk] = term_view[ii]
        aux_coeffs[kk] = merged_coeff
    end

    swapsums!(prop_cache)
    setactivesize!(prop_cache, lastactiveindex(prop_cache))

    return prop_cache
end
