
# =========================
# Vector propagation
# =========================

function PropagationBase.applymergetruncate!(gate::FermionicGate, prop_cache::MultiVectorMajoranaPropagationCache, theta; truncate_each_mr=nothing, kwargs...)
    # get the Majorana strings and coefficients corresponding to the fermionic gate
    #@show typeof(prop_cache)
    ms_rotations, coeffs, truncate_after_each_majrot = getmajoranarotations(gate, nsites(prop_cache))
    if !isnothing(truncate_each_mr)
        truncate_after_each_majrot = truncate_each_mr
    end

    #@show gate

    # iterate over individual Majorana rotations and apply them to the Majorana sum
    for (gate_ms, coeff) in zip(ms_rotations, coeffs)
        # check if gate is gaussian 
        is_gate_gaussian = _is_gaussian(gate_ms)


        pools = assign_pools(prop_cache)

        # multiply coefficient by 2 since `::MajoranaRotation` implements exp(-i * theta/2 * mstring)
        applytoall!(gate_ms, prop_cache, theta * coeff * 2.0; pools, kwargs...)

        # we need to merge 
        merge!(prop_cache; is_gaussian=is_gate_gaussian, kwargs...)

        # truncate after each Majorana rotation 
        if truncate_after_each_majrot
            truncate!(prop_cache; kwargs...)
        end
    end
    if !truncate_after_each_majrot
        truncate!(prop_cache; kwargs...)
    end

    return prop_cache
end


function PropagationBase.applytoall!(gate::MajoranaRotation{TT}, prop_cache::MultiVectorMajoranaPropagationCache, theta; pools, kwargs...) where {TT<:Integer,CT}
    #=for (weight_sector, vpropcache) in prop_caches(prop_cache)
        @show weight_sector, vpropcache.active_size
    end=#
    #=for (weight_sector, vpropcache) in prop_caches(prop_cache)
        applytoall!(gate, vpropcache, theta; is_gaussian=is_gaussian, _apply_function! = _applymajoranarotation_vm!, pool=pools[weight_sector], kwargs...)
        println("Applied MajoranaRotation to weight sector $weight_sector.")
        @show vpropcache.active_size
    end=#
    #@show length(prop_caches(prop_cache)[2])

    @sync for (weight_sector, vpropcache) in prop_caches(prop_cache)
        Threads.@spawn begin
            pool = pools[weight_sector]
            #println("Before applytoall ws=$weight_sector, activesize=$(activesize(vpropcache))")
            applytoall!(gate, vpropcache, theta; (_apply_function!)=_applymajoranarotation_vm!, pool=pool, kwargs...)
            #println("After applytoall ws=$weight_sector, activesize=$(activesize(vpropcache))")
        end
    end
    #pool = pools[first(keys(pools))]

    #=for (weight_sector, vpropcache) in prop_caches(prop_cache)
        applytoall!(gate, vpropcache, theta; (_apply_function!)=_applymajoranarotation_vm!, pool=nothing, kwargs...)
    end=#

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
    #tforeach(pool, eachindex(active_terms)) do ii
    for ii in eachindex(active_terms)
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
    #@show is_gaussian

    #=if is_gaussian
        println("Before merge ws=?, activesize=$(activesize(prop_cache))")
        merge!(prop_cache; kwargs...)
        println("After merge ws=?, activesize=$(activesize(prop_cache))")
    end=#

    return
end


function _is_gaussian(gate_ms::MajoranaRotation{TT}) where {TT<:Integer}
    return get_weight(gate_ms.ms_int) == 2
end

function assign_pools(prop_cache::MultiVectorMajoranaPropagationCache)
    max_threads = nthreads()
    pools = Dict{Int,ThreadPools.StaticPool}()
    sectors = collect(prop_caches(prop_cache))
    isempty(sectors) && return pools

    weights = Int[w for (w, _) in sectors]
    active_sizes = Float64[max(0, activesize(vpc)) for (_, vpc) in sectors]

    total_active = sum(active_sizes)
    nsec = length(sectors)

    # If all sectors are empty, give everyone one worker thread.
    if total_active == 0
        for w in weights
            pools[w] = ThreadPools.StaticPool(1:1)
        end
        return pools
    end

    # Valid, bounded allocation per sector: 1..max_threads.
    relative_sizes = active_sizes ./ total_active
    for i in eachindex(weights)
        n_chunks = clamp(round(Int, relative_sizes[i] * max_threads), 1, max_threads)
        pools[weights[i]] = ThreadPools.StaticPool(1:n_chunks)
    end

    return pools
end