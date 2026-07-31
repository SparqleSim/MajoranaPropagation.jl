
"""
    MajoranaRotation(ms::{TT}) where {TT<:Integer}
Basic structure to represent a Majorana rotation gate exp(-i * theta/2 * ms).
Defined by passing a Majorana string `ms` of even weight.
"""
struct MajoranaRotation{TT<:Integer} <: ParametrizedGate
    ms_int::TT
    function MajoranaRotation(ms::MajoranaString{TT}) where {TT<:Integer}
        # only even parity operations; the even weight is also a correctness invariant of the
        # specialized _commutes_evengate/_rotationproduct_evengate kernels used in applytoall!
        iseven(get_weight(ms)) || throw(ArgumentError("MajoranaRotation requires an even-weight Majorana string, got weight $(get_weight(ms))"))
        return new{TT}(ms.gammas)
    end
    function MajoranaRotation(ms_int::TT) where {TT<:Integer}
        iseven(get_weight(ms_int)) || throw(ArgumentError("MajoranaRotation requires an even-weight Majorana string, got weight $(get_weight(ms_int))"))
        return new{TT}(ms_int)
    end
end

"""
    FermionicGate(symbol::Symbol, sites::Vector{Int})
Structure to represent fermionic gates, constructed from a symbol. See `Constructors.jl` for supported symbols.
"""
struct FermionicRotation <: ParametrizedGate
    symbol::Symbol
    sites::Vector{Int}
end

function FermionicRotation(symbol::Symbol, site::Integer)
    return FermionicRotation(symbol, [site])
end

function FermionicRotation(symbol::Symbol, sites::Tuple)
    return FermionicRotation(symbol, collect(sites))
end


"""
    getmajoranarotations(gate::FermionicGate, n_sites::Integer)
Given a `FermionicGate`, returns the Majorana rotations and coefficients corresponding to it.
"""
function getmajoranarotations(gate::FermionicRotation, n_sites::Integer)
    # construct msum encoding the fermionic gate
    msum = MajoranaSum(n_sites, gate.symbol, gate.sites)
    TT = getinttype(nfermions(msum))
    truncate_after_each_majrot = !flag_non_number_preserving(gate.symbol)

    #remove coefficient associated to identity
    pop_id!(msum)

    rotations::Vector{MajoranaRotation{TT}} = []
    coefficients::Vector{Float64} = []
    for (ms, coeff) in msum
        push!(rotations, MajoranaRotation(ms))
        push!(coefficients, coeff)
    end

    return rotations, coefficients, truncate_after_each_majrot
end

function _applycos(coeff, cos_theta)
    return coeff * cos_theta
end
function _applysin(coeff, sin_theta)
    return coeff * sin_theta
end

"""
The splitting rule for exp(i theta gate_string / 2) ms exp(-i theta gate_string / 2) is
-) ms, if [gate_string, ms] = 0
-) cos(theta) ms - i sin(theta) ms * gate_string, if {gate_string, ms} = 0
"""
function PropagationBase.applytoall!(gate::MajoranaRotation, prop_cache::MajoranaPropagationCache, theta; kwargs...)
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    cos_val = cos(theta)
    sin_val = sin(theta)

    gate_int = gate.ms_int
    gate_int_ps = compute_parity_bits_and_shift(gate_int, 2 * nfermions(msum))
    omega_l_gate = omega_L_mult(gate_int)

    # loop over all Majorana strings and their coefficients in the Majorana sum
    for (ms_int, coeff) in msum
        if _commutes_evengate(ms_int, gate_int)
            # if the gate commutes with the Majorana string, do nothing
            continue
        end

        # else we know the gate will split the Majorana string into two
        coeff1 = _applycos(coeff, cos_val)
        new_ms, sign = _rotationproduct_evengate(ms_int, gate_int, gate_int_ps, omega_l_gate)
        coeff2 = _applysin(coeff, sin_val * sign)

        # set the coefficient of the original Majorana string
        set!(msum, ms_int, coeff1)

        # set the coefficient of the new Majorana string in the aux_psum
        # we can set the coefficient because MajoranaRotations create non-overlapping new Majorana strings
        set!(aux_msum, new_ms, coeff2)
    end

    return
end

function PropagationBase.applymergetruncate!(gate::FermionicRotation, prop_cache::AbstractMajoranaPropagationCache, theta; truncate_each_mr=nothing, kwargs...)
    # get the Majorana strings and coefficients corresponding to the fermionic gate
    ms_rotations, coeffs, truncate_after_each_majrot = getmajoranarotations(gate, nsites(prop_cache))
    if !isnothing(truncate_each_mr)
        truncate_after_each_majrot = truncate_each_mr
    end

    # iterate over individual Majorana rotations and apply them to the Majorana sum
    for (gate_ms, coeff) in zip(ms_rotations, coeffs)
        # multiply coefficient by 2 since `::MajoranaRotation` implements exp(-i * theta/2 * mstring)
        applytoall!(gate_ms, prop_cache, theta * coeff * 2.0; kwargs...)

        # merge the auxiliary Majorana sum into the original one and empty the auxiliary one
        # (for vector caches: the appended tail is merged into the sorted prefix tracked on
        # the sum; the gate string lets it sort the tail without a comparison sort)
        _mergeafterapply!(prop_cache, gate_ms.ms_int; kwargs...)

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

# bare MajoranaRotation gates take the same gate-aware merge as the rotations inside a
# FermionicRotation (the upstream generic applymergetruncate! would use the plain merge!,
# discarding the gate string that lets vector caches skip the comparison sort)
function PropagationBase.applymergetruncate!(gate::MajoranaRotation, prop_cache::AbstractMajoranaPropagationCache, theta; kwargs...)
    applytoall!(gate, prop_cache, theta; kwargs...)
    _mergeafterapply!(prop_cache, gate.ms_int; kwargs...)
    truncate!(prop_cache; kwargs...)
    return prop_cache
end


# ========== vector specializations ========== #

function PropagationBase.applytoall!(gate::MajoranaRotation, prop_cache::VectorMajoranaPropagationCache, theta; thread::Bool=true, kwargs...)

    if prop_cache.active_size == 0
        return prop_cache
    end

    n_old = prop_cache.active_size

    # get the Majorana string integer representation because the gate cannot be in the function when using GPU
    gate_ms = gate.ms_int
    gate_ms_ps = compute_parity_bits_and_shift(gate_ms, 2 * nfermions(prop_cache))
    omega_l_gate = omega_L_mult(gate_ms)

    # flag terms that anticommute with the gate
    anticommutesfunc(trm) = !_commutes_evengate(trm, gate_ms)
    flagterms!(anticommutesfunc, prop_cache; thread)

    # this runs a cumsum over the flags to get the indices
    flagstoindices!(prop_cache; thread)

    # the final index is the number of new terms
    n_noncommutes = lastactiveindex(prop_cache)

    # split off into the same array
    n_new = n_old + n_noncommutes

    # potential resize factor
    resize_factor = 2
    if capacity(prop_cache) < n_new
        resize!(prop_cache, n_new * resize_factor)
    end

    # does the branching logic
    _applymajoranarotation!(prop_cache, gate_ms, gate_ms_ps, omega_l_gate, theta; thread)

    # we now have n_new possibly duplicate Majorana strings in the array
    setactivesize!(prop_cache, n_new)

    return prop_cache
end

"""
The splitting rule for exp(i theta gate_string / 2) ms exp(-i theta gate_string / 2) is
-) ms, if [gate_string, ms] = 0
-) cos(theta) ms - i sin(theta) ms * gate_string, if {gate_string, ms} = 0
"""
function _applymajoranarotation!(prop_cache::VectorMajoranaPropagationCache, gate_ms::TT, gate_ms_ps::TT, omega_l_gate::Int, theta; thread::Bool=true) where {TT}

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
    AK.foreachindex(active_terms; max_tasks=maxtasks(thread), min_elems=_MIN_ELEMS_PER_TASK) do ii
        # here it anticommutes
        if flags[ii]
            term = terms[ii]
            coeff = coeffs[ii]

            coeff1 = coeff * cos_val
            new_term, sign = _rotationproduct_evengate(term, gate_ms, gate_ms_ps, omega_l_gate)
            coeff2 = coeff * sin_val * sign

            coeffs[ii] = coeff1

            terms[n+indices[ii]] = new_term
            coeffs[n+indices[ii]] = coeff2
        end
    end

    return
end