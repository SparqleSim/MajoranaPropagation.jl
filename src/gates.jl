
struct MajoranaRotation <: PauliPropagation.ParametrizedGate
    ms::MajoranaString
    function MajoranaRotation(ms::MajoranaString)
        @assert get_weight(ms) % 2 == 0 # only even parity operations
        return new(ms)
    end
end

function commutes(gate::MajoranaRotation, ms::MajoranaString)
    return commutes(gate.ms, ms)
end

struct FermionicGate <: PauliPropagation.ParametrizedGate
    symbol::Symbol
    sites::Vector{Int}
    nsites::Int
end

function FermionicGate(symbol::Symbol, site::Int, nfermions::Int)
    return FermionicGate(symbol, [site], nfermions)
end


function getmajoranarotations(gate::FermionicGate)
    # construct msum encoding the fermionic gate
    # multiply coefficient by 2 since exponential implements exp(-i * theta/2 * mstring)
    # TODO: fix where we multiply by 2
    msum = 2. * MajoranaSum(gate.nsites, gate.symbol, gate.sites)

    #remove coefficient associated to identity
    pop_id!(msum)

    return majoranas(msum), coefficients(msum)
end

function _applycos(coeff::CT, cos_theta::CT) where {CT}
    return coeff * cos_theta
end
function _applysin(coeff::CT, sin_theta::CT) where {CT}
    return coeff * sin_theta
end


function applytoall!(gate_int::TT, theta, msum::MajoranaSum{TT,CT}, aux_msum::MajoranaSum{TT,CT}; kwargs...) where {TT<:Integer,CT}
    cos_val = cos(theta)
    sin_val = sin(theta)

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

function applytoall!(gate::MajoranaRotation, theta, msum::MajoranaSum{TT,CT}, aux_msum::MajoranaSum{TT,CT}; kwargs...) where {TT<:Integer,CT}
    applytoall!(gate.ms.gammas, theta, msum, aux_msum; kwargs...)
end 

function applytoall!(gate::FermionicGate, theta, msum::MajoranaSum{TT,CT}, aux_msum::MajoranaSum{TT,CT}; kwargs...) where {TT<:Integer,CT}
    #get the Majorana strings and coefficients corresponding to the fermionic gate
    #@show gate.symbol gate.sites
    ms_rotations, coeffs = getmajoranarotations(gate)

    #iterate over individual Majorana rotations and apply them to the Majorana sum
    for (gate_ms, coeff) in zip(ms_rotations, coeffs)
        applytoall!(gate_ms, theta * coeff, msum, aux_msum; kwargs...)

        # merge the auxiliary Majorana sum into the original one and empty the auxiliary one
        msum, aux_msum = mergeandempty!(msum, aux_msum)

        # truncate after each Majorana rotation 
        checktruncationonall!(msum; kwargs...)
    end

    return 
end
