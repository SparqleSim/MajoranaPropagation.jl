using LinearAlgebra
using Bits

struct MajoranaString{TT<:Integer}
    nfermions::Int
    gammas::TT
end

function MajoranaString(nfermions::Int, indices::Vector{Int})
    TT = getinttype(nfermions)
    gammas::TT = 0
    for i in indices
        gammas = gammas | (typeof(gammas)(1) << (i - 1))
    end
    return MajoranaString(nfermions, gammas)
end

function MajoranaString(nfermions::Int, gammas::Int64)
    # Int64 is probably unwanted, lets make it the correct type
    TT = getinttype(nfermions)
    return MajoranaString(nfermions, convert(TT, gammas))
end

struct MajoranaSum{TT<:Integer,CT}
    nfermions::Int
    Majoranas::Dict{TT,CT}
end

function MajoranaSum(nfermions::Int, ::Type{CT}) where {CT}
    TT = getinttype(nfermions)
    return MajoranaSum(nfermions, Dict{TT,CT}())
end

function MajoranaSum(nfermions::Int, sumdict::Dict{Int,CT}) where {CT}
    TT = getinttype(nfermions)
    new_sumdict = Dict{TT,CT}()
    for (key, value) in sumdict
        new_sumdict[convert(TT, key)] = value
    end
    return MajoranaSum(nfermions, new_sumdict)
end

function MajoranaSum(nfermions::Int, ms::MajoranaString, value::CT) where {CT}
    TT = getinttype(nfermions)
    return MajoranaSum(nfermions, Dict{TT,CT}(ms.gammas => value))
end
function MajoranaSum(ms::MajoranaString, value::CT) where {CT}
    TT = getinttype(ms.nfermions)
    return MajoranaSum(ms.nfermions, Dict{TT,CT}(ms.gammas => value))
end

function MajoranaSum(nfermions::Int, ms_and_values::Vector{Tuple{CT,MajoranaString}}) where {CT}
    TT = getinttype(nfermions)
    sum_dict = Dict{TT,CT}()
    for (value, ms) in ms_and_values
        sum_dict[ms.gammas] = value
    end
    return MajoranaSum(nfermions, sum_dict)
end

function MajoranaSum(nfermions::Int, ms::MajoranaString, ::Type{CT}) where {CT}
    TT = getinttype(nfermions)
    return MajoranaSum(nfermions, Dict{TT,CT}(ms.gammas => CT(1.)))
end
function MajoranaSum(ms::MajoranaString, ::Type{CT}) where {CT}
    TT = getinttype(ms.nfermions)
    return MajoranaSum(ms.nfermions, Dict{TT,CT}(ms.gammas => CT(1.)))
end

function sum_add!(ms::MajoranaSum{TT,CT}, ms2::MajoranaString{TT}, value::CT) where {TT<:Integer,CT}
    sum_add!(ms, ms2.gammas, value)
end

function sum_add!(ms::MajoranaSum, ms2_gammas::TT, value::CT) where {TT<:Integer,CT}
    if haskey(ms.Majoranas, ms2_gammas)
        ms.Majoranas[ms2_gammas] += value
    else
        ms.Majoranas[ms2_gammas] = value
    end
end

function set!(ms::MajoranaSum{TT,CT}, ms2::MajoranaString, value::CT) where {TT<:Integer,CT}
    set!(ms, ms2.gammas, value)
    return
end

function set!(ms::MajoranaSum{TT,CT}, ms2::TT, value::CT) where {TT<:Integer,CT}
    ms.Majoranas[ms2] = value
    return
end

function get_coeffs(ms::MajoranaSum)
    return collect(values(ms.Majoranas))
end

function Base.delete!(ms::MajoranaSum, ms2::MajoranaString)
    delete!(ms.Majoranas, ms2.gammas)
end
function Base.delete!(ms::MajoranaSum, ms2_gammas::TT) where {TT<:Integer}
    delete!(ms.Majoranas, ms2_gammas)
end

function Base.pop!(ms::MajoranaSum, ms2_gammas::TT) where {TT<:Integer}
    return pop!(ms.Majoranas, ms2_gammas, 0.)
end

function Base.length(ms::MajoranaSum)
    return length(ms.Majoranas)
end

function Base.mergewith!(merge, msum1::MajoranaSum, msum2::MajoranaSum)
    mergewith!(merge, msum1.Majoranas, msum2.Majoranas)
    return msum1
end

function Base.empty!(msum::MajoranaSum)
    empty!(msum.Majoranas)
    return msum
end

function Base.show(io::IO, ms::MajoranaString)
    print(io, "$(reverse(string(ms.gammas; base=2, pad=2 * ms.nfermions)))")
end

function Base.show(io::IO, ms::MajoranaSum)
    max_display = 8
    print(io, "MajoranaSum with $(length(ms)) term(s):(")
    for (i, (mstring, coeff)) in enumerate(ms.Majoranas)
        if i <= max_display
            print(io, "\n")
            print(io, "    $(coeff) * $(reverse(string(mstring; base=2, pad=2 * ms.nfermions)))")
        else
            print(io, "\n    ...")
            break
        end
    end
    print(io, ")")
end


function coefftype(msum::MajoranaSum{TT,CT}) where {TT,CT}
    return CT
end

function similar(msum::MajoranaSum)
    return MajoranaSum(msum.nfermions, coefftype(msum))
end

Base.iterate(msum::MajoranaSum, state=1) = iterate(msum.Majoranas, state)

function get_weight(ms::MajoranaString)
    return get_weight(ms.gammas)
end
function get_weight(gammas::TT) where {TT<:Integer}
    return Bits.weight(gammas)
end

function int_to_binary_array(var::TT, N::Int) where {TT<:Integer}
    var_string = string(var; base=2, pad=N)

    bits_array = zeros(Int, N)
    for i in 1:N
        bits_array[i] = parse(Int, var_string[N+1-i])
    end
    return bits_array
end

function compute_parity_bits_and_shift(u::TT, Nbits::Int) where {TT<:Integer}
    p::TT = 0
    parity::TT = 0
    for k in TT(0):TT(Nbits - 2)
        if (u >> k) & 1 == 1
            parity ⊻= 1
        end
        p |= parity << (k + 1)
    end
    return p
end

function omega_L_mult(ms1::MajoranaString, ms2::MajoranaString)
    return omega_L_mult(ms1.gammas, ms2.gammas, 2 * ms1.nfermions)
end

function omega_L_mult(ms1::TT, ms2::TT, Nbits) where {TT<:Integer}
    return mod(Bits.weight(ms1 & compute_parity_bits_and_shift(ms2, Nbits)), 2)
end

function omega_L_mult(ms::TT) where {TT<:Integer}
    wms = get_weight(ms)
    return mod((wms^2 - wms) / 2, 2)
end

function omega_L_mult(ms::MajoranaString)
    return omega_L_mult(ms.gammas)
end

function omega_mult(ms1::MajoranaString, ms2::MajoranaString)
    return omega_mult(ms1.gammas, ms2.gammas)
end

function omega_mult(gammas1::TT, gammas2::TT) where {TT<:Integer}
    w1 = get_weight(gammas1)
    w2 = get_weight(gammas2)
    return mod(w1 * w2 - get_weight(gammas1 & gammas2), 2)
end

function omega_mult(ms::MajoranaString)
    return omega_L_mult(ms, ms)
end

function Base.:(==)(ms1::MajoranaString, ms2::MajoranaString)
    if ms1.nfermions != ms2.nfermions
        return false
    end
    return ms1.gammas == ms2.gammas
end

function Base.:(==)(ms1::MajoranaSum, ms2::MajoranaSum)
    if ms1.nfermions != ms2.nfermions
        return false
    end
    return ms1.Majoranas == ms2.Majoranas
end

function mstring_additon(ms1::TT, ms2::TT) where {TT<:Integer}
    return ms1 ⊻ ms2
end
function Base.:(+)(ms1::MajoranaString, ms2::MajoranaString)
    if ms1.nfermions != ms2.nfermions
        throw(ArgumentError("Majorana strings must have the same length, but have lengths $(ms1.nfermions) and $(ms2.nfermions)"))
    end
    return MajoranaString(ms1.nfermions, mstring_additon(ms1.gammas, ms2.gammas))
end

function Base.:(*)(msum1::MajoranaSum, msum2::MajoranaSum)
    res = MajoranaSum(msum1.nfermions, coefftype(msum1))
    for (ms1, coeff1) in msum1.Majoranas
        for (ms2, coeff2) in msum2.Majoranas
            prefactor, ms3 = ms_mult(MajoranaString(msum1.nfermions, ms1), MajoranaString(msum2.nfermions, ms2))
            @assert imag(prefactor) ≈ 0
            prefactor = real(prefactor)
            sum_add!(res, ms3, prefactor * tonumber(coeff1) * tonumber(coeff2))
        end
    end
    return res
end

function fprefactor(g1::TT, g2::TT) where {TT<:Integer}
    return omega_L_mult(g1) * omega_L_mult(g2) + omega_mult(g1, g2) * (omega_L_mult(g1) + omega_L_mult(g2) + 1)
end

function fprefactor(ms1::MajoranaString, ms2::MajoranaString)
    return fprefactor(ms1.gammas, ms2.gammas)
end

function ms_mult(ms1::MajoranaString, ms2::MajoranaString)
    if ms1.nfermions != ms2.nfermions
        throw(ArgumentError("Majorana strings must have the same length, but have lengths $(ms1.nfermions) and $(ms2.nfermions)"))
    end
    prefactor, result = ms_mult(ms1.gammas, ms2.gammas, 2 * ms1.nfermions)
    return prefactor, MajoranaString(ms1.nfermions, result)
end

function ms_mult(ms1::TT, ms2::TT, n_fermions::Integer) where {TT<:Integer}
    result = mstring_additon(ms1, ms2) # result = ms1 + ms2
    prefactor = (-1)^(omega_L_mult(ms1, ms2, 2 * n_fermions) + fprefactor(ms1, ms2))
    if mod(omega_mult(ms1, ms2), 2) == 1
        return 1im * prefactor, result
    end
    return prefactor, result
end

function commutes(ms1::MajoranaString, ms2::MajoranaString)
    return commutes(ms1.gammas, ms2.gammas)
end

function commutes(gammas1::Integer, gammas2::Integer)
    return mod(omega_mult(gammas1, gammas2), 2) == 0
end


function norm(msum::MajoranaSum, L=2)
    if length(msum) == 0
        return 0.0
    end
    return LinearAlgebra.norm((coeff for coeff in get_coeffs(msum)), L)
end

function commutator(msum1::MajoranaSum, msum2::MajoranaSum)
    res = MajoranaSum(msum1.nfermions, typeof(1.1im))
    for (ms1, coeff1) in msum1.Majoranas
        for (ms2, coeff2) in msum2.Majoranas
            if commutes(ms1, ms2)
                continue
            end
            prefactor, ms3 = ms_mult(MajoranaString(msum1.nfermions, ms1), MajoranaString(msum2.nfermions, ms2))
            sum_add!(res, ms3, prefactor * tonumber(coeff1) * tonumber(coeff2))
        end
    end
    return res
end

function pop_id!(msum::MajoranaSum)
    if haskey(msum.Majoranas, 0)
        delete!(msum.Majoranas, 0)
    end
    return
end

function fock_filter(msum::MajoranaSum)
    clean_res = similar(msum)
    singles_filter = create_max_single_filter(msum.nfermions)
    for (ms, coeff) in msum.Majoranas
        if compute_max_single(ms, 2 * msum.nfermions, singles_filter) > 0
            continue
        end
        set!(clean_res, ms, coeff)
    end
    return clean_res
end

function overlap_with_fock(msum::MajoranaSum, fock_state; add_pref=0.)
    res = 0.
    non_contributing = 0
    singles_filter = create_max_single_filter(msum.nfermions)
    for (ms, coeff) in msum.Majoranas
        res += fockevaluate(ms, coeff, singles_filter, fock_state)
    end
    return res + add_pref
end

function fockevaluate(ms::TT, coeff, singles_filter, fock_state) where {TT<:Integer}
    if compute_max_single(ms, 2, singles_filter) > 0
        return 0.
    end
    num_pref = 0
    for site in fock_state
        num_pref += (ms >> (2 * site - 1)) & 1
    end
    ms_w = get_weight(ms)
    sign = (1im)^omega_L_mult(ms) * (1im)^(ms_w / 2) * (-1)^num_pref
    return tonumber(coeff) * sign
end

function overlap_with_fock_spinful(mslist, up_sites_with_particle, down_sites_with_particle, nsites; add_pref=0.)
    fock_state = []
    for up_site in up_sites_with_particle
        push!(fock_state, 2 * up_site - 1)
    end
    for down_site in down_sites_with_particle
        push!(fock_state, 2 * down_site)
    end
    return overlap_with_fock(mslist, fock_state; add_pref=add_pref)
end
