#using PauliPropagation.PropagationBase
#import PauliPropagation.PropagationBase: mainsum, auxsum
#using MajoranaPropagation: AbstractMajoranaPropagationCache, VectorMajoranaPropagationCache, AbstractMajoranaSum, majoranas
#import AcceleratedKernels
#const AK = AcceleratedKernels

struct ImaginaryMajoranaRotation{TT<:Integer} <: ParametrizedGate
    ms_int::TT
    function ImaginaryMajoranaRotation(ms::MajoranaString{TT}) where {TT<:Integer}
        iseven(get_weight(ms)) || throw(ArgumentError("ImaginaryMajoranaRotation requires an even-weight Majorana string, got weight $(get_weight(ms))"))
        return new{TT}(ms.gammas)
    end
    function ImaginaryMajoranaRotation(ms_int::TT) where {TT<:Integer}
        iseven(get_weight(ms_int)) || throw(ArgumentError("ImaginaryMajoranaRotation requires an even-weight Majorana string, got weight $(get_weight(ms_int))"))
        return new{TT}(ms_int)
    end
end


struct ImaginaryFermionicRotation <: ParametrizedGate
    symbol::Symbol
    sites::Vector{Int}
end

function ImaginaryFermionicRotation(symbol::Symbol, site::Integer)
    return ImaginaryFermionicRotation(symbol, [site])
end

function PauliPropagation._toheisenberg(gate::Union{ImaginaryFermionicRotation,ImaginaryMajoranaRotation}, τ)
    throw(error("$(typeof(gate)) gates are currently not defined in the Heisenberg picture."))
end

function PauliPropagation._toschrodinger(gate::Union{ImaginaryFermionicRotation,ImaginaryMajoranaRotation}, τ)
    return gate, τ
end

"""
    getmajoranarotations(gate::ImaginaryFermionicGate, n_sites::Integer)
Given a `ImaginaryFermionicGate`, returns the Majorana rotations and coefficients corresponding to it.
"""
function getmajoranarotations(gate::ImaginaryFermionicRotation, n_sites::Integer)
    # construct msum encoding the fermionic gate
    msum = MajoranaSum(n_sites, gate.symbol, gate.sites)
    TT = getinttype(nfermions(msum))
    truncate_after_each_majrot = !flag_non_number_preserving(gate.symbol)

    #remove coefficient associated to identity
    pop_id!(msum)

    rotations::Vector{ImaginaryMajoranaRotation{TT}} = []
    coefficients::Vector{Float64} = []
    for (ms, coeff) in msum
        push!(rotations, ImaginaryMajoranaRotation(ms))
        push!(coefficients, coeff)
    end

    return rotations, coefficients, truncate_after_each_majrot
end

"""
Implement exp(- beta fermionic_gate / 2) msum exp(- beta fermionic_gate / 2) for imaginary time evolution
"""
function PropagationBase.applymergetruncate!(gate::ImaginaryFermionicRotation, prop_cache::AbstractMajoranaPropagationCache, beta; truncate_each_mr=nothing, normalize_coeffs=true, kwargs...)
    # get the Majorana strings and coefficients corresponding to the fermionic gate
    ms_rotations, coeffs, truncate_after_each_majrot = getmajoranarotations(gate, nsites(prop_cache))
    if !isnothing(truncate_each_mr)
        truncate_after_each_majrot = truncate_each_mr
    end

    # iterate over individual Majorana rotations and apply them to the Majorana sum
    for (gate_ms, coeff) in zip(ms_rotations, coeffs)
        applytoall!(gate_ms, prop_cache, beta * coeff; kwargs...)

        # merge the auxiliary Majorana sum into the original one and empty the auxiliary one
        # (for vector caches: the appended tail is gate ⊻ (ascending commuting terms), so the
        # same XOR-sorted tail merge as in the real-time path applies)
        _mergeafterapply!(prop_cache, gate_ms.ms_int; kwargs...)

        # normalize coefficients to preserve state normalization
        if normalize_coeffs
            mult!(prop_cache, 1 / getmergedcoeff(mainsum(prop_cache), 0))
        end

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

"""
Implement exp(- beta majorana_rotation / 2) msum exp(- beta majorana_rotation / 2) for imaginary time evolution
For a Majorana string ms, the splitting rule for exp(- beta majorana_rotation / 2) ms exp(- beta majorana_rotation / 2) is
-) ms, if {majorana_rotation, ms} = 0
-) cosh(beta) ms - sinh(beta)  majorana_rotation * ms, if [majorana_rotation, ms] = 0
"""
function PropagationBase.applytoall!(gate::ImaginaryMajoranaRotation, prop_cache::MajoranaPropagationCache, beta; kwargs...)
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    cosh_val = cosh(beta)
    sinh_val = -sinh(beta)

    gate_int = gate.ms_int
    gate_int_ps = compute_parity_bits_and_shift(gate_int, 2 * nfermions(msum))
    omega_l_gate = omega_L_mult(gate_int)

    # loop over all Majorana strings and their coefficients in the Majorana sum
    for (ms_int, coeff) in msum
        if _commutes_evengate(ms_int, gate_int)

            # the imaginary gate will split the Majorana string into two
            coeff1 = coeff * cosh_val
            new_ms, sign = _rotationproduct_evengate_commuting(ms_int, gate_int, gate_int_ps, omega_l_gate)
            coeff2 = coeff * sinh_val * sign

            # set the coefficient of the original Majorana string
            set!(msum, ms_int, coeff1)

            # set the coefficient of the new Majorana string in the aux_psum
            set!(aux_msum, new_ms, coeff2)
        else
            continue
        end
    end

    return

end

# ========== vector specializations ========== #
"""
Implement exp(- beta majorana_rotation / 2) msum exp(- beta majorana_rotation / 2) for imaginary time evolution
"""
function PropagationBase.applytoall!(gate::ImaginaryMajoranaRotation, prop_cache::VectorMajoranaPropagationCache, beta; thread::Bool=true, kwargs...)

    if prop_cache.active_size == 0
        return prop_cache
    end

    n_old = prop_cache.active_size

    # get the Majorana string integer representation because the gate cannot be in the function when using GPU
    gate_ms = gate.ms_int
    # gate invariants hoisted once per gate; gate_ms has even weight by construction
    gate_ms_ps = compute_parity_bits_and_shift(gate_ms, 2 * nfermions(prop_cache))
    omega_l_gate = omega_L_mult(gate_ms)

    # in imaginary time we split upon commutation
    commutesfunc(trm) = _commutes_evengate(trm, gate_ms)
    PropagationBase.flagterms!(commutesfunc, prop_cache; thread)

    # this runs a cumsum over the flags to get the indices
    PropagationBase.flagstoindices!(prop_cache; thread)

    # the final index is the number of new terms
    n_commutes = PropagationBase.lastactiveindex(prop_cache)

    # split off into the same array
    n_new = n_old + n_commutes

    # potential resize factor
    resize_factor = 1.5
    if capacity(prop_cache) < n_new
        resize!(prop_cache, round(Int, n_new * resize_factor))
    end

    # does the branching logic
    _applyimaginarymajoranarotation!(prop_cache, gate_ms, gate_ms_ps, omega_l_gate, beta; thread)

    # we now have n_new possibly duplicate Majorana strings in the array
    PropagationBase.setactivesize!(prop_cache, n_new)

    return prop_cache
end

"""
For a Majorana string ms, the splitting rule for exp(- beta majorana_rotation / 2) ms exp(- beta majorana_rotation / 2) is
-) ms, if {majorana_rotation, ms} = 0
-) cosh(beta) ms - sinh(beta)  majorana_rotation * ms, if [majorana_rotation, ms] = 0
"""
function _applyimaginarymajoranarotation!(prop_cache::VectorMajoranaPropagationCache, gate_ms::TT, gate_ms_ps::TT, omega_l_gate::Int, beta; thread::Bool=true) where {TT}

    # pre-compute the sine and cosine values because they are used for every Majorana string that does not commute with the gate
    cosh_val = cosh(beta)
    sinh_val = -sinh(beta)

    n = PropagationBase.activesize(prop_cache)
    n_max = n + PropagationBase.lastactiveindex(prop_cache)

    active_terms = PropagationBase.activeterms(prop_cache)

    # full-length terms so we can write new terms at the end
    terms = majoranas(PropagationBase.mainsum(prop_cache))
    coeffs = coefficients(PropagationBase.mainsum(prop_cache))
    @assert length(terms) >= n_max "VectorMajoranaPropagationCache terms array is not large enough to hold new terms."
    @assert length(coeffs) >= n_max "VectorMajoranaPropagationCache coeffs array is not large enough to hold new coeffs."

    flags = PropagationBase.activeflags(prop_cache)
    indices = PropagationBase.activeindices(prop_cache)

    # branching pattern for Majorana rotations
    AK.foreachindex(active_terms; max_tasks=_maxtasks(thread)) do ii
        # here it commutes
        if flags[ii]
            term = terms[ii]
            coeff = coeffs[ii]

            coeff1 = coeff * cosh_val
            new_term, sign = _rotationproduct_evengate_commuting(term, gate_ms, gate_ms_ps, omega_l_gate)
            coeff2 = coeff * sinh_val * sign

            coeffs[ii] = coeff1

            terms[n+indices[ii]] = new_term
            coeffs[n+indices[ii]] = coeff2
        end
    end

    return
end

