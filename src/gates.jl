
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
    # mstrings given in the "Schroedinger order"
    mstrings::Vector{MajoranaString}
    relative_signs::Vector{Float64}

    function FermionicGate(mstrings::Vector{MajoranaString}, relative_signs::Vector{Float64})
        @assert length(mstrings) == length(relative_signs)
        return new(mstrings, relative_signs)
    end
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

    gate_int::TT = gate.ms.gammas
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

function applytoall!(gate::FermionicGate, theta, msum::MajoranaSum{TT,CT}, aux_msum::MajoranaSum{TT,CT}; kwargs...) where {TT<:Integer,CT}
    cos_val = cos(theta)

    for (gate_ms, rel_sign) in zip(reverse(gate.mstrings), reverse(gate.relative_signs))
        gate_int::TT = gate_ms.gammas
        sin_val = sin(theta * rel_sign)
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
        #after each gate, merge the aux_msum back into the main Majorana sum
        msum, aux_msum = mergeandempty!(msum, aux_msum)
    end
    return
end
