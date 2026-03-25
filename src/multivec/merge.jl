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
    _flagstoindices_pool!(activeindices(prop_cache), activeflags(prop_cache); pool)
    _mergegroups_pool!(prop_cache; pool)
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

function _flagstoindices_pool!(dst_indices, flags; pool)
    n = length(flags)
    if n == 0
        return dst_indices
    end

    nblocks = min(n, 8 * Threads.nthreads())
    block_size = cld(n, nblocks)
    block_sums = zeros(Int, nblocks)
    block_offsets = zeros(Int, nblocks)

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
    tforeach(pool, 1:nblocks) do b
        lo = (b - 1) * block_size + 1
        hi = min(b * block_size, n)
        local_count = block_offsets[b]
        for ii in lo:hi
            local_count += flags[ii] ? 1 : 0
            dst_indices[ii] = local_count
        end
    end

    return dst_indices
end

function _mergegroups_pool!(prop_cache::VectorMajoranaPropagationCache; pool)
    term_view = activeterms(prop_cache)
    coeffs = activecoeffs(prop_cache)
    aux_terms = activeauxterms(prop_cache)
    aux_coeffs = activeauxcoeffs(prop_cache)
    flags = activeflags(prop_cache)
    indices = activeindices(prop_cache)
    active_size = activesize(prop_cache)

    group_starts = Int[]
    for ii in eachindex(flags)
        if flags[ii]
            push!(group_starts, ii)
        end
    end

    tforeach(pool, eachindex(group_starts)) do kk
        ii = group_starts[kk]
        end_idx = kk < length(group_starts) ? group_starts[kk + 1] - 1 : active_size

        CT = typeof(coeffs[ii])
        merged_coeff = zero(CT)
        for jj in ii:end_idx
            merged_coeff = mergefunc(merged_coeff, coeffs[jj])
        end

        aux_terms[indices[ii]] = term_view[ii]
        aux_coeffs[indices[ii]] = merged_coeff
    end

    swapsums!(prop_cache)
    setactivesize!(prop_cache, lastactiveindex(prop_cache))

    return prop_cache
end
