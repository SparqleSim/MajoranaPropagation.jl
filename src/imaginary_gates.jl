using PauliPropagation.PropagationBase
import PauliPropagation.PropagationBase: mainsum, auxsum
using MajoranaPropagation: AbstractMajoranaPropagationCache, VectorMajoranaPropagationCache, AbstractMajoranaSum, majoranas
import AcceleratedKernels
const AK = AcceleratedKernels

struct ImaginaryMajoranaRotation{TT<:Integer} <: ParametrizedGate
    ms::MajoranaString{TT}
    function ImaginaryMajoranaRotation(ms::MajoranaString{TT}) where {TT<:Integer}
        @assert get_weight(ms) % 2 == 0 # only even parity operations
        return new{TT}(ms)
    end
end


struct ImaginaryFermionicGate <: ParametrizedGate
    symbol::Symbol
    sites::Vector{Int}
end

function ImaginaryFermionicGate(symbol::Symbol, site::Integer)
    return ImaginaryFermionicGate(symbol, [site])
end


MajoranaPropagation.nfermions(prop_cache::AbstractMajoranaPropagationCache) = nfermions(mainsum(prop_cache))

function PropagationBase.applymergetruncate!(gate::ImaginaryFermionicGate, prop_cache::AbstractMajoranaPropagationCache, theta; kwargs...)
    # get the Majorana strings and coefficients corresponding to the fermionic gate
    ms_rotations, coeffs = MajoranaPropagation.getmajoranarotations(gate, nsites(prop_cache))

    # iterate over individual Majorana rotations and apply them to the Majorana sum
    for (gate_ms, coeff) in zip(ms_rotations, coeffs)
        # multiply coefficient by 2 since `::MajoranaRotation` implements exp(-i * theta/2 * mstring)
        MajoranaPropagation.applytoall!(gate_ms, prop_cache, theta * coeff * 2.0; kwargs...)

        # merge the auxiliary Majorana sum into the original one and empty the auxiliary one
        MajoranaPropagation.merge!(prop_cache; kwargs...)

        # truncate after each Majorana rotation 
        MajoranaPropagation.truncate!(prop_cache; kwargs...)
    end

    return prop_cache
end

function MajoranaPropagation.getmajoranarotations(gate::ImaginaryFermionicGate, n_sites::Integer)
    # construct msum encoding the fermionic gate
    msum = MajoranaSum(n_sites, gate.symbol, gate.sites)

    #remove coefficient associated to identity
    MajoranaPropagation.pop_id!(msum)

    rotations::Vector{ImaginaryMajoranaRotation} = []
    coefficients::Vector{Float64} = []
    for (ms, coeff) in msum
        push!(rotations, ImaginaryMajoranaRotation(MajoranaString(nfermions(msum), ms)))
        push!(coefficients, coeff)
    end

    return rotations, coefficients
end

# ========== vector specializations ========== #

function MajoranaPropagation.applytoall!(gate::ImaginaryMajoranaRotation, prop_cache::VectorMajoranaPropagationCache, theta; kwargs...)

    if prop_cache.active_size == 0
        return prop_cache
    end

    n_old = prop_cache.active_size

    # get the Majorana string integer representation because the gate cannot be in the function when using GPU
    gate_ms = gate.ms.gammas

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

