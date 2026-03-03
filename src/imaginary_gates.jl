#using PauliPropagation.PropagationBase
#import PauliPropagation.PropagationBase: mainsum, auxsum
#using MajoranaPropagation: AbstractMajoranaPropagationCache, VectorMajoranaPropagationCache, AbstractMajoranaSum, majoranas
#import AcceleratedKernels
#const AK = AcceleratedKernels

struct ImaginaryMajoranaRotation{TT<:Integer} <: ParametrizedGate
    ms_int::TT
    function ImaginaryMajoranaRotation(ms::MajoranaString{TT}) where {TT<:Integer}
        @assert get_weight(ms) % 2 == 0 # only even parity operations
        return new{TT}(ms.gammas)
    end
    function ImaginaryMajoranaRotation(ms_int::TT) where {TT<:Integer}
        @assert get_weight(ms_int) % 2 == 0 # only even parity operations
        return new{TT}(ms_int)
    end
end


struct ImaginaryFermionicGate <: ParametrizedGate
    symbol::Symbol
    sites::Vector{Int}
end

function ImaginaryFermionicGate(symbol::Symbol, site::Integer)
    return ImaginaryFermionicGate(symbol, [site])
end

function PauliPropagation._toheisenberg(gate::Union{ImaginaryFermionicGate,ImaginaryMajoranaRotation}, τ)
    throw(error("$(typeof(gate)) gates are currently not defined in the Heisenberg picture."))
end

function PauliPropagation._toschrodinger(gate::Union{ImaginaryFermionicGate,ImaginaryMajoranaRotation}, τ)
    return gate, τ
end

"""
    getmajoranarotations(gate::ImaginaryFermionicGate, n_sites::Integer)
Given a `ImaginaryFermionicGate`, returns the Majorana rotations and coefficients corresponding to it.
"""
function getmajoranarotations(gate::ImaginaryFermionicGate, n_sites::Integer)
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


function PropagationBase.applymergetruncate!(gate::ImaginaryFermionicGate, prop_cache::AbstractMajoranaPropagationCache, theta; truncate_each_mr=nothing, normalize_coeffs=true, kwargs...)
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
        merge!(prop_cache; kwargs...)

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

function PropagationBase.applytoall!(gate::ImaginaryMajoranaRotation, prop_cache::MajoranaPropagationCache, theta; kwargs...)
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    cosh_val = cosh(theta)
    sinh_val = sinh(theta)

    gate_int = gate.ms_int

    # loop over all Majorana strings and their coefficients in the Majorana sum
    for (ms_int, coeff) in msum
        if commutes(gate_int, ms_int)

            # the imaginary gate will split the Majorana string into two
            coeff1 = coeff * cosh_val
            sign, new_ms = ms_mult(gate_int, ms_int, nfermions(msum))
            coeff2 = coeff * sinh_val * real(sign) # TODO: there might be a -1 missing

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

function PropagationBase.applytoall!(gate::ImaginaryMajoranaRotation, prop_cache::VectorMajoranaPropagationCache, theta; kwargs...)

    if prop_cache.active_size == 0
        return prop_cache
    end

    n_old = prop_cache.active_size

    # get the Majorana string integer representation because the gate cannot be in the function when using GPU
    gate_ms = gate.ms_int

    # in imaginary time we split upon commutation
    commutesfunc(trm) = MajoranaPropagation.commutes(trm, gate_ms)
    PropagationBase.flagterms!(commutesfunc, prop_cache)

    # this runs a cumsum over the flags to get the indices
    PropagationBase.flagstoindices!(prop_cache)

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
    _applyimaginarymajoranarotation!(prop_cache, gate_ms, theta)

    # we now have n_new possibly duplicate Majorana strings in the array
    PropagationBase.setactivesize!(prop_cache, n_new)

    return prop_cache
end

function _applyimaginarymajoranarotation!(prop_cache::VectorMajoranaPropagationCache, gate_ms::TT, theta) where {TT}

    # pre-compute the sine and cosine values because they are used for every Majorana string that does not commute with the gate
    cosh_val = cosh(theta)
    sinh_val = sinh(theta)

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
    AK.foreachindex(active_terms) do ii
        # here it anticommutes
        if flags[ii]
            term = terms[ii]
            coeff = coeffs[ii]

            coeff1 = coeff * cosh_val
            sign, new_term = ms_mult(gate_ms, term, nfermions(prop_cache))
            coeff2 = coeff * sinh_val * real(sign) # TODO: there might be a -1 missing

            coeffs[ii] = coeff1

            terms[n+indices[ii]] = new_term
            coeffs[n+indices[ii]] = coeff2
        end
    end

    return
end

