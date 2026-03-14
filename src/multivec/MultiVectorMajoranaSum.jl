using Base.Threads

"""
    MultiVectorMajoranaSum

A weighted-sector Majorana sum where each weight sector is stored as a
`VectorMajoranaSum`. Internally this is a dictionary from Majorana weight
to a `VectorMajoranaSum` restricted to that weight.

This combines:

- the contiguous vector storage of `VectorMajoranaSum`, and
- the parallel weight-sector merging strategy of `MajoranaSumMulti`.
"""
struct MultiVectorMajoranaSum{TT<:Integer,CT} <: AbstractMajoranaSum
    nsites::Int
    is_spinful::Bool
    MultiMajoranas::Dict{Int,VectorMajoranaSum{Vector{TT},Vector{CT}}}
end

"""Construct a `MultiVectorMajoranaSum` from a dense `MajoranaSum`."""
function MultiVectorMajoranaSum(msum::MajoranaSum{TT,CT}) where {TT<:Integer,CT}
    nsites = msum.nsites
    is_spinful = msum.is_spinful

    multimajs = Dict{Int,VectorMajoranaSum{Vector{TT},Vector{CT}}}()

    for (ms_int, coeff) in msum.Majoranas
        weight = get_weight(ms_int)
        vms = get!(multimajs, weight) do
            VectorMajoranaSum(nsites, is_spinful, TT[], CT[])
        end
        push!(vms.terms, ms_int)
        push!(vms.coeffs, coeff)
    end

    return MultiVectorMajoranaSum{TT,CT}(nsites, is_spinful, multimajs)
end

PropagationBase.storage(msum::MultiVectorMajoranaSum) = msum.MultiMajoranas

Base.keys(msum::MultiVectorMajoranaSum) = Base.keys(msum.MultiMajoranas)

function Base.length(msum::MultiVectorMajoranaSum{TT,CT}) where {TT<:Integer,CT}
    total_strings = 0
    for (_w, vms) in msum.MultiMajoranas
        total_strings += length(vms)
    end
    return total_strings
end

function coefftype(::MultiVectorMajoranaSum{TT,CT}) where {TT<:Integer,CT}
    return CT
end

PropagationBase.nsites(msum::MultiVectorMajoranaSum) = msum.nsites

function nfermions(msum::MultiVectorMajoranaSum)
    if msum.is_spinful
        return 2 * msum.nsites
    else
        return msum.nsites
    end
end

"""
    similar(msum::MultiVectorMajoranaSum, W::Int)

Create an empty `MultiVectorMajoranaSum` which has support only on the
sectors with weights `W-2`, `W`, and `W+2`. This mirrors the behaviour
of `similar(::MajoranaSumMulti, W)` and is used by the propagation cache.
"""
function similar(msum::MultiVectorMajoranaSum{TT,CT}, W::Int) where {TT<:Integer,CT}
    multimajs = Dict{Int,VectorMajoranaSum{Vector{TT},Vector{CT}}}()
    for w in (W - 2, W, W + 2)
        multimajs[w] = VectorMajoranaSum(msum.nsites, msum.is_spinful, TT[], CT[])
    end
    return MultiVectorMajoranaSum{TT,CT}(msum.nsites, msum.is_spinful, multimajs)
end

"""
    show_stats(msum::MultiVectorMajoranaSum)

Print statistics on the number of strings per weight sector.
"""
function show_stats(msum::MultiVectorMajoranaSum{TT,CT}) where {TT<:Integer,CT}
    total_strings = length(msum)
    if total_strings == 0
        println("MultiVectorMajoranaSum is empty.")
        return
    end
    sorted_keys = sort(collect(keys(msum.MultiMajoranas)))
    for weight_key in sorted_keys
        nstrings = length(msum.MultiMajoranas[weight_key])
        println("Weight $weight_key: $nstrings strings ($(round(100.0 * nstrings / total_strings))%)")
    end
    println("Total strings: $total_strings")
end


# =========================
# Propagation cache
# =========================

mutable struct MultiVectorMajoranaPropagationCache{MVMS<:MultiVectorMajoranaSum} <: AbstractMajoranaPropagationCache
    main_msum::MVMS
    aux_msum::Dict{Integer,MVMS}
end

PropagationBase.mainsum(prop_cache::MultiVectorMajoranaPropagationCache) = prop_cache.main_msum
PropagationBase.auxsum(prop_cache::MultiVectorMajoranaPropagationCache) = prop_cache.aux_msum
nfermions(prop_cache::MultiVectorMajoranaPropagationCache) = nfermions(mainsum(prop_cache))

function PropagationBase.setmainsum!(prop_cache::MultiVectorMajoranaPropagationCache, msum::MultiVectorMajoranaSum)
    prop_cache.main_msum = msum
    return prop_cache
end

function PropagationBase.setauxsum!(
    prop_cache::MultiVectorMajoranaPropagationCache,
    aux_msum::Dict{Integer,<:MultiVectorMajoranaSum},
)
    prop_cache.aux_msum = aux_msum
    return prop_cache
end

"""
    MultiVectorMajoranaPropagationCache(msum::MultiVectorMajoranaSum)

Create a propagation cache for a `MultiVectorMajoranaSum`. For each
weight sector in the main sum we allocate an auxiliary multi-vector sum
with support in the neighbouring weight sectors.
"""
function MultiVectorMajoranaPropagationCache(multimsum::MultiVectorMajoranaSum{TT,CT}) where {TT<:Integer,CT}
    all_aux_msums::Dict{Integer,MultiVectorMajoranaSum{TT,CT}} = Dict()
    for weight_key in keys(multimsum)
        all_aux_msums[weight_key] = similar(multimsum, weight_key)
    end
    return MultiVectorMajoranaPropagationCache(multimsum, all_aux_msums)
end

PropagationBase.PropagationCache(multimsum::MultiVectorMajoranaSum) = MultiVectorMajoranaPropagationCache(multimsum)


# =========================
# Vector-style rotation per weight sector
# =========================

"""
    _push_new_term!(mvsum, ms_int, coeff)

Helper to push a new Majorana string and coefficient into the
appropriate weight sector of a `MultiVectorMajoranaSum`.
"""
function _push_new_term!(mvsum::MultiVectorMajoranaSum{TT,CT}, ms_int::TT, coeff::CT) where {TT<:Integer,CT}
    weight = get_weight(ms_int)
    vms = get!(mvsum.MultiMajoranas, weight) do
        VectorMajoranaSum(mvsum.nsites, mvsum.is_spinful, TT[], CT[])
    end
    push!(vms.terms, ms_int)
    push!(vms.coeffs, coeff)
end

function _append_sector_terms!(
    mvsum::MultiVectorMajoranaSum{TT,CT},
    weight::Int,
    terms::Vector{TT},
    coeffs::Vector{CT},
) where {TT<:Integer,CT}
    if isempty(terms)
        return
    end
    vms = get!(mvsum.MultiMajoranas, weight) do
        VectorMajoranaSum(mvsum.nsites, mvsum.is_spinful, TT[], CT[])
    end
    append!(vms.terms, terms)
    append!(vms.coeffs, coeffs)
    return
end

@inline function _chunk_bounds(n::Int, nchunks::Int, chunk::Int)
    start_idx = fld((chunk - 1) * n, nchunks) + 1
    end_idx = fld(chunk * n, nchunks)
    return start_idx, end_idx
end

mutable struct SectorThreadPool{TT<:Integer,CT}
    ranges::Vector{UnitRange{Int}}
    local_terms_m2::Vector{Vector{TT}}
    local_coeffs_m2::Vector{Vector{CT}}
    local_terms_w::Vector{Vector{TT}}
    local_coeffs_w::Vector{Vector{CT}}
    local_terms_p2::Vector{Vector{TT}}
    local_coeffs_p2::Vector{Vector{CT}}
end

function _build_sector_thread_pool(::Type{TT}, ::Type{CT}, n::Int, n_chunks::Int) where {TT<:Integer,CT}
    nch = min(max(n_chunks, 1), max(n, 1))
    ranges = UnitRange{Int}[]
    for chunk in 1:nch
        s, e = _chunk_bounds(n, nch, chunk)
        if s <= e
            push!(ranges, s:e)
        end
    end
    nworkers = length(ranges)
    return SectorThreadPool{TT,CT}(
        ranges,
        [TT[] for _ in 1:nworkers],
        [CT[] for _ in 1:nworkers],
        [TT[] for _ in 1:nworkers],
        [CT[] for _ in 1:nworkers],
        [TT[] for _ in 1:nworkers],
        [CT[] for _ in 1:nworkers],
    )
end

function _allocate_chunks_by_size(
    sector_sizes::Vector{Int},
    max_chunks::Int;
    min_chunk_size::Int=2048,
)
    nsectors = length(sector_sizes)
    if nsectors == 0
        return Int[]
    end

    chunks = zeros(Int, nsectors)

    # Small sectors are executed sequentially on one shared worker.
    small = findall(<(min_chunk_size), sector_sizes)
    large = setdiff(collect(1:nsectors), small)

    if isempty(large)
        return chunks
    end

    # Reserve one worker for the sequential small-sector lane when needed.
    reserved_for_small = isempty(small) ? 0 : 1
    budget = max(1, max_chunks - reserved_for_small)

    large_sizes = [sector_sizes[i] for i in large]
    total_large = sum(large_sizes)
    if total_large == 0
        chunks[large[1]] = 1
        return chunks
    end

    # Proportional allocation by sector-size fraction.
    raw = [budget * (s / total_large) for s in large_sizes]
    alloc = floor.(Int, raw)
    rem = raw .- alloc

    # Largest remainder fill.
    leftover = budget - sum(alloc)
    if leftover > 0
        order = sortperm(rem; rev=true)
        for j in 1:leftover
            alloc[order[(j - 1) % length(order) + 1]] += 1
        end
    end

    # Ensure each large sector gets at least one chunk.
    missing = findall(==(0), alloc)
    if !isempty(missing)
        donors = sortperm(alloc; rev=true)
        for m in missing
            for d in donors
                if alloc[d] > 1
                    alloc[d] -= 1
                    alloc[m] += 1
                    break
                end
            end
        end
    end

    for (k, idx) in enumerate(large)
        chunks[idx] = min(alloc[k], sector_sizes[idx])
    end

    return chunks
end

"""
    _applymajoranarotation!(vms, aux_msum, gate_int, theta, n_fermions; merge_sector=false, weight_key)

Apply a single `MajoranaRotation` (specified by its integer representation `gate_int`)
to all strings inside a single weight sector `vms`. New strings are written into
the auxiliary multi-vector sum `aux_msum`.
"""
function _applymajoranarotation!(
    vms::VectorMajoranaSum{Vector{TT},Vector{CT}},
    aux_msum::MultiVectorMajoranaSum{TT,CT},
    gate_int,
    theta,
    n_fermions;
    merge_sector::Bool=false,
    weight_key::Int,
    thread_pool::Union{Nothing,SectorThreadPool{TT,CT}}=nothing,
    kwargs...,
) where {TT<:Integer,CT}
    cos_val = cos(theta)
    sin_val = sin(theta)

    terms, coeffs = PropagationBase.storage(vms)

    if isnothing(thread_pool)
        @inbounds for i in eachindex(terms, coeffs)
            term = terms[i]
            coeff = coeffs[i]
            if commutes(gate_int, term)
                continue
            end
            coeff1 = coeff * cos_val
            sign, new_term = ms_mult(gate_int, term, n_fermions)
            coeff2 = coeff * sin_val * real((-1im) * sign)
            coeffs[i] = coeff1
            _push_new_term!(aux_msum, new_term, coeff2)
        end
    else
        ranges = thread_pool.ranges
        local_terms_m2 = thread_pool.local_terms_m2
        local_coeffs_m2 = thread_pool.local_coeffs_m2
        local_terms_w = thread_pool.local_terms_w
        local_coeffs_w = thread_pool.local_coeffs_w
        local_terms_p2 = thread_pool.local_terms_p2
        local_coeffs_p2 = thread_pool.local_coeffs_p2

        for chunk_id in eachindex(ranges)
            empty!(local_terms_m2[chunk_id])
            empty!(local_coeffs_m2[chunk_id])
            empty!(local_terms_w[chunk_id])
            empty!(local_coeffs_w[chunk_id])
            empty!(local_terms_p2[chunk_id])
            empty!(local_coeffs_p2[chunk_id])
        end

        @sync for chunk_id in eachindex(ranges)
            @spawn begin
                range = ranges[chunk_id]

                chunk_terms_m2 = local_terms_m2[chunk_id]
                chunk_coeffs_m2 = local_coeffs_m2[chunk_id]
                chunk_terms_w = local_terms_w[chunk_id]
                chunk_coeffs_w = local_coeffs_w[chunk_id]
                chunk_terms_p2 = local_terms_p2[chunk_id]
                chunk_coeffs_p2 = local_coeffs_p2[chunk_id]

                @inbounds for i in range
                    term = terms[i]
                    coeff = coeffs[i]
                    if commutes(gate_int, term)
                        continue
                    end
                    coeff1 = coeff * cos_val
                    sign, new_term = ms_mult(gate_int, term, n_fermions)
                    coeff2 = coeff * sin_val * real((-1im) * sign)
                    coeffs[i] = coeff1
                    new_weight = get_weight(new_term)
                    if new_weight == weight_key - 2
                        push!(chunk_terms_m2, new_term)
                        push!(chunk_coeffs_m2, coeff2)
                    elseif new_weight == weight_key + 2
                        push!(chunk_terms_p2, new_term)
                        push!(chunk_coeffs_p2, coeff2)
                    else
                        push!(chunk_terms_w, new_term)
                        push!(chunk_coeffs_w, coeff2)
                    end
                end
            end
        end

        for chunk_id in eachindex(ranges)
            _append_sector_terms!(aux_msum, weight_key - 2, local_terms_m2[chunk_id], local_coeffs_m2[chunk_id])
            _append_sector_terms!(aux_msum, weight_key, local_terms_w[chunk_id], local_coeffs_w[chunk_id])
            _append_sector_terms!(aux_msum, weight_key + 2, local_terms_p2[chunk_id], local_coeffs_p2[chunk_id])
        end
    end

    if merge_sector
        # Merge back only the same-weight sector
        if haskey(aux_msum.MultiMajoranas, weight_key)
            aux_sector = aux_msum.MultiMajoranas[weight_key]
            vms = _merge_vectorsums!(vms, aux_sector)
            # empty the auxiliary sector
            empty!(aux_sector.terms)
            empty!(aux_sector.coeffs)
        end
    end

    return
end

# =========================
# Merging logic
# =========================

"""
    _merge_vectorsums!(dest, src)

Merge `src` into `dest`, adding coefficients of identical strings.
The `src` vectors are emptied.
"""
function _merge_vectorsums!(
    dest::VectorMajoranaSum{Vector{TT},Vector{CT}},
    src::VectorMajoranaSum{Vector{TT},Vector{CT}},
) where {TT<:Integer,CT}
    d_terms, d_coeffs = PropagationBase.storage(dest)
    s_terms, s_coeffs = PropagationBase.storage(src)

    nd = length(d_terms)
    ns = length(s_terms)

    if ns == 0
        return dest
    elseif nd == 0
        # move src into dest: copy into existing storage
        append!(d_terms, s_terms)
        append!(d_coeffs, s_coeffs)
        empty!(s_terms)
        empty!(s_coeffs)
        return dest
    end

    # sort both sector vectors by term using AcceleratedKernels.sort!
    sort!(dest)
    sort!(src)

    total = nd + ns
    merged_terms = Vector{TT}(undef, total)
    merged_coeffs = Vector{CT}(undef, total)

    i = 1
    j = 1
    k = 1

    while i <= nd && j <= ns
        td = d_terms[i]
        ts = s_terms[j]
        if td == ts
            merged_terms[k] = td
            merged_coeffs[k] = d_coeffs[i] + s_coeffs[j]
            i += 1
            j += 1
        elseif td < ts
            merged_terms[k] = td
            merged_coeffs[k] = d_coeffs[i]
            i += 1
        else
            merged_terms[k] = ts
            merged_coeffs[k] = s_coeffs[j]
            j += 1
        end
        k += 1
    end

    while i <= nd
        merged_terms[k] = d_terms[i]
        merged_coeffs[k] = d_coeffs[i]
        i += 1
        k += 1
    end

    while j <= ns
        merged_terms[k] = s_terms[j]
        merged_coeffs[k] = s_coeffs[j]
        j += 1
        k += 1
    end

    final_len = k - 1
    resize!(merged_terms, final_len)
    resize!(merged_coeffs, final_len)

    empty!(d_terms)
    empty!(d_coeffs)
    append!(d_terms, merged_terms)
    append!(d_coeffs, merged_coeffs)

    empty!(s_terms)
    empty!(s_coeffs)

    return dest
end

function Base.merge!(prop_cache::MultiVectorMajoranaPropagationCache; merge_sector::Bool=true, kwargs...)
    prop_cache = _merge!(prop_cache; merge_sector, kwargs...)
    return prop_cache
end

function _merge!(
    prop_cache::MultiVectorMajoranaPropagationCache{MVMS};
    merge_sector::Bool=true,
    kwargs...,
) where {TT<:Integer,CT,MVMS<:MultiVectorMajoranaSum{TT,CT}}
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    sorted_keys = sort(collect(keys(msum)))

    # ΔW = 0 merges
    if merge_sector
        @threads for weight_key in sorted_keys
            if haskey(aux_msum, weight_key)
                main_dict = msum.MultiMajoranas
                aux_dict = aux_msum[weight_key].MultiMajoranas
                if haskey(aux_dict, weight_key)
                    if haskey(main_dict, weight_key)
                        _merge_vectorsums!(main_dict[weight_key], aux_dict[weight_key])
                    else
                        main_dict[weight_key] = aux_dict[weight_key]
                    end
                    empty!(aux_dict[weight_key].terms)
                    empty!(aux_dict[weight_key].coeffs)
                end
            end
        end
    end

    # ΔW = +2 merges
    @threads for weight_key in sorted_keys
        if haskey(aux_msum, weight_key)
            aux_dict = aux_msum[weight_key].MultiMajoranas
            target_key = weight_key + 2
            if haskey(aux_dict, target_key) && !isempty(aux_dict[target_key].terms)
                main_dict = msum.MultiMajoranas
                if haskey(main_dict, target_key)
                    _merge_vectorsums!(main_dict[target_key], aux_dict[target_key])
                    empty!(aux_dict[target_key].terms)
                    empty!(aux_dict[target_key].coeffs)
                else
                    # move sector into main and replace aux sector with an empty one
                    main_dict[target_key] = aux_dict[target_key]
                    aux_dict[target_key] = VectorMajoranaSum(msum.nsites, msum.is_spinful, eltype(aux_dict[target_key].terms)[], eltype(aux_dict[target_key].coeffs)[])
                end
            end
        end
    end

    # ΔW = -2 merges
    @threads for weight_key in sorted_keys
        if haskey(aux_msum, weight_key)
            aux_dict = aux_msum[weight_key].MultiMajoranas
            target_key = weight_key - 2
            if haskey(aux_dict, target_key) && !isempty(aux_dict[target_key].terms)
                main_dict = msum.MultiMajoranas
                if haskey(main_dict, target_key)
                    _merge_vectorsums!(main_dict[target_key], aux_dict[target_key])
                    empty!(aux_dict[target_key].terms)
                    empty!(aux_dict[target_key].coeffs)
                else
                    # move sector into main and replace aux sector with an empty one
                    main_dict[target_key] = aux_dict[target_key]
                    aux_dict[target_key] = VectorMajoranaSum(msum.nsites, msum.is_spinful, eltype(aux_dict[target_key].terms)[], eltype(aux_dict[target_key].coeffs)[])
                end
            end
        end
    end

    # remove empty sectors
    for weight_key in sort(collect(keys(msum.MultiMajoranas)))
        if isempty(msum.MultiMajoranas[weight_key].terms)
            delete!(msum.MultiMajoranas, weight_key)
        end
    end

    setmainsum!(prop_cache, msum)
    setauxsum!(prop_cache, aux_msum)

    return prop_cache
end


# =========================
# applytoall! specialization
# =========================

function PropagationBase.applytoall!(
    gate::MajoranaRotation,
    prop_cache::MultiVectorMajoranaPropagationCache,
    theta;
    merge_sector::Bool=false,
    kwargs...,
)
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    gate_int = gate.ms_int
    msum_keys = collect(keys(msum))
    if isempty(msum_keys)
        return prop_cache
    end

    sector_sizes = [length(msum.MultiMajoranas[w]) for w in msum_keys]
    n_chunks_per_sector = _allocate_chunks_by_size(sector_sizes, nthreads(); min_chunk_size=2048)
    small_sector_idxs = findall(==(0), n_chunks_per_sector)
    large_sector_idxs = findall(>(0), n_chunks_per_sector)

    # Ensure aux entries exist before starting parallel tasks.
    for weight_key in msum_keys
        if !haskey(aux_msum, weight_key)
            aux_msum[weight_key] = similar(msum, weight_key)
        end
    end

    sector_pools = Dict{Int,Any}()
    for iw in large_sector_idxs
        weight_key = msum_keys[iw]
        vms = msum.MultiMajoranas[weight_key]
        terms, coeffs = PropagationBase.storage(vms)
        TT = eltype(terms)
        CT = eltype(coeffs)
        sector_pools[weight_key] = _build_sector_thread_pool(TT, CT, length(terms), n_chunks_per_sector[iw])
    end

    @sync begin
        # One shared worker for all small sectors, processed sequentially.
        if !isempty(small_sector_idxs)
            @spawn begin
                for iw in small_sector_idxs
                    weight_key = msum_keys[iw]
                    vms = msum.MultiMajoranas[weight_key]
                    aux_entry = aux_msum[weight_key]
                    _applymajoranarotation!(
                        vms,
                        aux_entry,
                        gate_int,
                        theta,
                        nfermions(msum);
                        merge_sector=merge_sector,
                        weight_key=weight_key,
                        thread_pool=nothing,
                        kwargs...,
                    )
                end
            end
        end

        # Large sectors run concurrently, each with proportional chunk count.
        for iw in large_sector_idxs
            @spawn begin
                weight_key = msum_keys[iw]
                vms = msum.MultiMajoranas[weight_key]
                aux_entry = aux_msum[weight_key]
                _applymajoranarotation!(
                    vms,
                    aux_entry,
                    gate_int,
                    theta,
                    nfermions(msum);
                    merge_sector=merge_sector,
                    weight_key=weight_key,
                    thread_pool=sector_pools[weight_key],
                    kwargs...,
                )
            end
        end
    end

    return prop_cache
end


# =========================
# Truncation
# =========================

function PropagationBase.truncate!(
    prop_cache::MultiVectorMajoranaPropagationCache;
    max_weight::Real=Inf,
    min_abs_coeff=1e-10,
    max_unpaired::Real=Inf,
    max_freq::Real=Inf,
    max_sins::Real=Inf,
    unpaired_mask=nothing,
    customtruncfunc=nothing,
    kwargs...,
)
    if isnothing(unpaired_mask)
        unpaired_mask = create_unpaired_mask(nfermions(mainsum(prop_cache)))
    end

    function truncfunc(mstr, coeff)
        is_truncated = false
        if PauliPropagation.truncatemincoeff(coeff, min_abs_coeff)
            is_truncated = true
        elseif truncateunpaired(mstr, max_unpaired, unpaired_mask)
            is_truncated = true
        elseif truncatemajoranaweight(mstr, max_weight)
            is_truncated = true
        elseif PauliPropagation.truncatefrequency(coeff, max_freq)
            is_truncated = true
        elseif PauliPropagation.truncatesins(coeff, max_sins)
            is_truncated = true
        elseif !isnothing(customtruncfunc) && customtruncfunc(mstr, coeff)
            is_truncated = true
        end
        return is_truncated
    end

    msum = mainsum(prop_cache)

    # loop over weight sectors and compact in-place
    for weight_key in collect(keys(msum.MultiMajoranas))
        vms = msum.MultiMajoranas[weight_key]
        terms, coeffs = PropagationBase.storage(vms)

        write_idx = 0
        for i in eachindex(terms, coeffs)
            ms_int = terms[i]
            coeff = coeffs[i]
            if truncfunc(ms_int, coeff)
                continue
            end
            write_idx += 1
            terms[write_idx] = ms_int
            coeffs[write_idx] = coeff
        end

        if write_idx == 0
            delete!(msum.MultiMajoranas, weight_key)
        else
            resize!(terms, write_idx)
            resize!(coeffs, write_idx)
        end
    end

    setmainsum!(prop_cache, msum)

    return
end

