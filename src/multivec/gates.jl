using TimerOutputs
# =========================
# Vector propagation
# =========================

function PropagationBase.applymergetruncate!(gate::FermionicGate, prop_cache::MultiVectorMajoranaPropagationCache, theta; truncate_each_mr=nothing, to::TimerOutput, kwargs...)
    # get the Majorana strings and coefficients corresponding to the fermionic gate
    ms_rotations, coeffs, truncate_after_each_majrot = getmajoranarotations(gate, nsites(prop_cache))
    if !isnothing(truncate_each_mr)
        truncate_after_each_majrot = truncate_each_mr
    end

    # iterate over individual Majorana rotations and apply them to the Majorana sum
    for (gate_ms, coeff) in zip(ms_rotations, coeffs)
        # check if gate is gaussian 
        is_gate_gaussian = _is_gaussian(gate_ms)

        pools = assign_pools(prop_cache)

        # multiply coefficient by 2 since `::MajoranaRotation` implements exp(-i * theta/2 * mstring)
        @timeit to "applytoall!" applytoall!(gate_ms, prop_cache, theta * coeff * 2.0; pools, kwargs...)

        # we need to merge 
        @timeit to "merge!" merge!(prop_cache; is_gaussian=is_gate_gaussian, to, pools, kwargs...)

        # truncate after each Majorana rotation 
        if truncate_after_each_majrot
            @timeit to "truncate!" truncate!(prop_cache; kwargs...)
        end
    end
    if !truncate_after_each_majrot
        @timeit to "truncate!" truncate!(prop_cache; kwargs...)
    end

    return prop_cache
end


function PropagationBase.applytoall!(gate::MajoranaRotation{TT}, prop_cache::MultiVectorMajoranaPropagationCache, theta; pools, kwargs...) where {TT<:Integer,CT}
    @sync for (weight_sector, vpropcache) in prop_caches(prop_cache)
        Threads.@spawn begin
            pool = pools[weight_sector]
            applytoall!(gate, vpropcache, theta; (_apply_function!)=_applymajoranarotation_vm!, pool=pool, kwargs...)
        end
    end

    return prop_cache
end

function _applymajoranarotation_vm!(prop_cache::VectorMajoranaPropagationCache, gate_ms::TT, theta; pool, kwargs...) where {TT}
    # pre-compute the sine and cosine values because they are used for every Majorana string that does not commute with the gate
    cos_val = cos(theta)
    sin_val = sin(theta)
    n_fermions = nfermions(prop_cache)

    n = activesize(prop_cache)
    n_max = n + lastactiveindex(prop_cache)

    active_terms = activeterms(prop_cache)

    # full-length terms so we can write new terms at the end
    terms = majoranas(mainsum(prop_cache))
    coeffs = coefficients(mainsum(prop_cache))
    @assert length(terms) >= n_max "VectorMajoranaPropagationCache terms array is not large enough to hold new terms."
    @assert length(coeffs) >= n_max "VectorMajoranaPropagationCache coeffs array is not large enough to hold new coeffs."

    flags = activeflags(prop_cache)
    indices = activeindices(prop_cache)

    # branching pattern for Majorana rotations
    tforeach(pool, eachindex(active_terms)) do ii
        # here it anticommutes
        if flags[ii]
            term = terms[ii]
            coeff = coeffs[ii]

            coeff1 = coeff * cos_val
            sign, new_term = ms_mult(gate_ms, term, n_fermions)
            coeff2 = coeff * sin_val * real((-1im) * sign)

            coeffs[ii] = coeff1

            terms[n+indices[ii]] = new_term
            coeffs[n+indices[ii]] = coeff2
        end
    end

    return
end


function _is_gaussian(gate_ms::MajoranaRotation{TT}) where {TT<:Integer}
    return get_weight(gate_ms.ms_int) == 2
end

function assign_pools(prop_cache::MultiVectorMajoranaPropagationCache; max_threads=nthreads())
    pools = Dict{Int,ThreadPools.StaticPool}()
    sectors = collect(prop_caches(prop_cache))
    isempty(sectors) && return pools

    weights = Int[w for (w, _) in sectors]

    if max_threads == 1
        for w in weights
            pools[w] = ThreadPools.StaticPool(1:1)
        end
        return pools
    end

    active_sizes = Float64[max(0, activesize(vpc)) for (_, vpc) in sectors]
    total_active = sum(active_sizes)

    # If all sectors are empty, give everyone one worker thread.
    if total_active == 0
        for w in weights
            pools[w] = ThreadPools.StaticPool(1:1)
        end
        return pools
    end

    return assign_pools(Dict(zip(weights, active_sizes ./ total_active)); max_threads=max_threads)
end

function assign_pools(active_sizes_percent::Dict{Int,Float64}; max_threads=nthreads())
    pools = Dict{Int,ThreadPools.StaticPool}()

    weights = collect(keys(active_sizes_percent))
    relative_sizes = collect(values(active_sizes_percent))

    sorting_indices = sortperm(relative_sizes, rev=true)
    weights = weights[sorting_indices]
    relative_sizes = relative_sizes[sorting_indices]

    prev_thread = 2

    for i in eachindex(weights)
        n_chunks = max(1, floor(Int, relative_sizes[i] * max_threads))
        if prev_thread > max_threads
            pools[weights[i]] = ThreadPools.StaticPool(1:1)
        else
            if prev_thread + n_chunks > max_threads
                pools[weights[i]] = ThreadPools.StaticPool(prev_thread:max_threads)
                prev_thread += n_chunks
            else
                pools[weights[i]] = ThreadPools.StaticPool(prev_thread:(prev_thread+n_chunks-1))
                prev_thread += n_chunks
            end
        end
    end

    return pools
end