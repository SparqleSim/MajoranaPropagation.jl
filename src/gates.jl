
"""
    MajoranaRotation(ms::MajoranaString{TT}) where {TT<:Integer}
Basic structure to represent a Majorana rotation gate exp(-i * theta/2 * ms).
Defined by passing a Majorana string `ms` of even weight.
"""
struct MajoranaRotation{TT<:Integer} <: ParametrizedGate
    ms::MajoranaString{TT}
    function MajoranaRotation(ms::MajoranaString{TT}) where {TT<:Integer}
        @assert get_weight(ms) % 2 == 0 # only even parity operations
        return new{TT}(ms)
    end
end

"""
    FermionicGate(symbol::Symbol, sites::Vector{Int})
Structure to represent fermionic gates, constructed from a symbol. See `Constructors.jl` for supported symbols.
"""
struct FermionicGate <: ParametrizedGate
    symbol::Symbol
    sites::Vector{Int}
end

function FermionicGate(symbol::Symbol, site::Integer)
    return FermionicGate(symbol, [site])
end


"""
    getmajoranarotations(gate::FermionicGate, n_sites::Integer)
Given a `FermionicGate`, returns the Majorana rotations and coefficients corresponding to it.
"""
function getmajoranarotations(gate::FermionicGate, n_sites::Integer)
    # construct msum encoding the fermionic gate
    msum = MajoranaSum(n_sites, gate.symbol, gate.sites)

    #remove coefficient associated to identity
    pop_id!(msum)

    rotations::Vector{MajoranaRotation} = []
    coefficients::Vector{Float64} = []
    for (ms, coeff) in msum
        push!(rotations, MajoranaRotation(MajoranaString(nfermions(msum), ms)))
        push!(coefficients, coeff)
    end

    return rotations, coefficients
end

function _applycos(coeff, cos_theta)
    return coeff * cos_theta
end
function _applysin(coeff, sin_theta)
    return coeff * sin_theta
end


function PropagationBase.applytoall!(gate::MajoranaRotation, prop_cache::MajoranaPropagationCache, theta; kwargs...)
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    cos_val = cos(theta)
    sin_val = sin(theta)

    gate_int = gate.ms.gammas

    # loop over all Majorana strings and their coefficients in the Majorana sum
    for (ms_int, coeff) in msum
        if commutes(gate_int, ms_int)
            # if the gate commutes with the Majorana string, do nothing
            continue
        end

        # else we know the gate will split the Majorana string into two
        coeff1 = _applycos(coeff, cos_val)
        sign, new_ms = ms_mult(gate_int, ms_int, nfermions(msum))
        coeff2 = _applysin(coeff, sin_val * real((-1im) * sign))

        # set the coefficient of the original Majorana string
        set!(msum, ms_int, coeff1)

        # set the coefficient of the new Majorana string in the aux_psum
        # we can set the coefficient because MajoranaRotations create non-overlapping new Majorana strings
        set!(aux_msum, new_ms, coeff2)
    end

    return
end

function PropagationBase.applymergetruncate!(gate::FermionicGate, prop_cache::MajoranaPropagationCache, theta; kwargs...)
    # get the Majorana strings and coefficients corresponding to the fermionic gate
    ms_rotations, coeffs = getmajoranarotations(gate, nsites(prop_cache))

    # iterate over individual Majorana rotations and apply them to the Majorana sum
    for (gate_ms, coeff) in zip(ms_rotations, coeffs)
        # multiply coefficient by 2 since exponential implements exp(-i * theta/2 * mstring)
        applytoall!(gate_ms, prop_cache, theta * coeff * 2.0; kwargs...)

        # merge the auxiliary Majorana sum into the original one and empty the auxiliary one
        merge!(prop_cache; kwargs...)

        # truncate after each Majorana rotation 
        truncate!(prop_cache; kwargs...)
    end

    return prop_cache
end
