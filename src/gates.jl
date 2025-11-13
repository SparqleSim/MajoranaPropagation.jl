
"""
    MajoranaRotation(ms::MajoranaString{TT}) where {TT<:Integer}
Basic structure to represent a Majorana rotation gate exp(-i * theta/2 * ms).
Defined by passing a Majorana string `ms` of even weight.
"""
struct MajoranaRotation{TT<:Integer} <: PauliPropagation.ParametrizedGate
    ms::MajoranaString{TT}
    function MajoranaRotation(ms::MajoranaString{TT}) where {TT<:Integer}
        @assert get_weight(ms) % 2 == 0 # only even parity operations
        return new{TT}(ms)
    end
end

"""
    FermionicGate(symbol::Symbol, sites::Vector{Int}, nfermions::Int)
Structure to represent fermionic gates, constructed from a symbol. See `Constructors.jl` for supported symbols.
"""
struct FermionicGate <: PauliPropagation.ParametrizedGate
    symbol::Symbol
    sites::Vector{Int}
    nsites::Int
end

function FermionicGate(symbol::Symbol, site::Int, nfermions::Int)
    return FermionicGate(symbol, [site], nfermions)
end


"""
    getmajoranarotations(gate::FermionicGate)
Given a `FermionicGate`, returns the Majorana rotations and coefficients corresponding to it.
"""
function getmajoranarotations(gate::FermionicGate)
    # construct msum encoding the fermionic gate
    msum = MajoranaSum(gate.nsites, gate.symbol, gate.sites)

    #remove coefficient associated to identity
    pop_id!(msum)

    rotations::Vector{MajoranaRotation} = []
    coefficients::Vector{Float64} = []
    for (ms, coeff) in msum
        push!(rotations, MajoranaRotation(MajoranaString(msum.nfermions, ms)))
        push!(coefficients, coeff)
    end

    return rotations, coefficients
end

function _applycos(coeff::CT, cos_theta::CT) where {CT}
    return coeff * cos_theta
end
function _applysin(coeff::CT, sin_theta::CT) where {CT}
    return coeff * sin_theta
end


function applytoall!(gate::MajoranaRotation, theta, msum::MajoranaSum{TT,CT}, aux_msum::MajoranaSum{TT,CT}; kwargs...) where {TT<:Integer,CT}
    cos_val = cos(theta)
    sin_val = sin(theta)

    gate_int = gate.ms.gammas

    # loop over all Majorana strings and their coefficients in the Majorana sum
    for (ms_int, coeff) in msum.Majoranas
        if commutes(gate_int, ms_int)
            # if the gate commutes with the Majorana string, do nothing
            continue
        end

        # else we know the gate will split the Majorana string into two
        coeff1 = _applycos(coeff, cos_val)
        sign, new_ms = ms_mult(gate_int, ms_int, msum.nfermions)
        coeff2 = _applysin(coeff, sin_val * real((-1im) * sign))

        # set the coefficient of the original Majorana string
        set!(msum, ms_int, coeff1)

        # set the coefficient of the new Majorana string in the aux_psum
        # we can set the coefficient because MajoranaRotations create non-overlapping new Majorana strings
        set!(aux_msum, new_ms, coeff2)
    end

    return
end

function applymergetruncate!(gate::FermionicGate, msum::MajoranaSum{TT,CT}, aux_msum::MajoranaSum{TT,CT}, thetas, param_idx; kwargs...) where {TT<:Integer,CT}
    # get the Majorana strings and coefficients corresponding to the fermionic gate
    ms_rotations, coeffs = getmajoranarotations(gate)

    # get the current parameter
    theta = thetas[param_idx]

    # iterate over individual Majorana rotations and apply them to the Majorana sum
    for (gate_ms, coeff) in zip(ms_rotations, coeffs)
        # multiply coefficient by 2 since exponential implements exp(-i * theta/2 * mstring)
        applytoall!(gate_ms, theta * coeff * 2.0, msum, aux_msum; kwargs...)

        # merge the auxiliary Majorana sum into the original one and empty the auxiliary one
        msum, aux_msum = mergeandempty!(msum, aux_msum)

        # truncate after each Majorana rotation 
        checktruncationonall!(msum; kwargs...)
    end

    # decrement parameter index during back propagation
    return msum, aux_msum, param_idx - 1
end
