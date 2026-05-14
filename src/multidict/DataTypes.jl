import Base: keys
import PauliPropagation.PropagationBase: add!

struct MajoranaSumMulti{TT<:Integer,CT} <: AbstractMajoranaSum
    nsites::Int
    is_spinful::Bool
    MultiMajoranas::Vector{Dict{TT,CT}}
    level_mapper::Function
    n_levels::Int
end

function MajoranaSumMulti(msum::MajoranaSum{TT,CT}, level_mapper::Function, n_levels::Int) where {TT<:Integer,CT}
    nsites = msum.nsites
    is_spinful = msum.is_spinful
    multimajs = [Dict{TT,CT}() for _ in 1:n_levels]

    for (ms_int, coeff) in msum.Majoranas
        ms_level = level_mapper(ms_int)
        multimajs[ms_level][ms_int] = coeff
    end
    return MajoranaSumMulti{TT,CT}(nsites, is_spinful, multimajs, level_mapper, n_levels)
end

# default level mapper
function MajoranaSumMulti(msum::MajoranaSum{TT,CT}) where {TT<:Integer,CT}
    n_levels = nthreads()
    level_mapper = ms -> rem(rem(ms * 11, n_levels) + countweight(ms), n_levels) + 1
    return MajoranaSumMulti(msum, level_mapper, n_levels)
end

function Base.keys(msum::MajoranaSumMulti)
    return Base.keys(msum.MultiMajoranas)
end

# set ms of certain weight assuming dict weight is already present in msum
function set!(msum::MajoranaSumMulti{TT,CT}, level_key, ms_int::TT, coeff::CT) where {TT<:Integer,CT}
    msum.MultiMajoranas[level_key][ms_int] = coeff
end

function add!(msum::MajoranaSumMulti{TT,CT}, level_key, ms_int::TT, coeff::CT) where {TT<:Integer,CT}
    dict = msum.MultiMajoranas[level_key]
    dict[ms_int] = get(dict, ms_int, zero(CT)) + coeff
    return msum
end

function similar(msum::MajoranaSumMulti{TT,CT}) where {TT<:Integer,CT}
    out_vec = [Dict{TT,CT}() for _ in 1:length(msum.MultiMajoranas)]
    return MajoranaSumMulti(msum.nsites, msum.is_spinful, out_vec, msum.level_mapper, msum.n_levels)
end

function PropagationBase.terms(msum::MajoranaSumMulti{TT,CT}) where {TT<:Integer,CT}
    all_ms = Vector{TT}(undef, length(msum))
    idx = 1
    for wi in eachindex(msum.MultiMajoranas)
        dict = msum.MultiMajoranas[wi]
        for ms_int in keys(dict)
            all_ms[idx] = ms_int
            idx += 1
        end
    end
    return all_ms
end

function PropagationBase.coefficients(msum::MajoranaSumMulti{TT,CT}) where {TT<:Integer,CT}
    all_coeffs = Vector{CT}(undef, length(msum))
    idx = 1
    for wi in eachindex(msum.MultiMajoranas)
        dict = msum.MultiMajoranas[wi]
        for coeff in values(dict)
            all_coeffs[idx] = coeff
            idx += 1
        end
    end
    return all_coeffs
end

function coefftype(msum::MajoranaSumMulti{TT,CT}) where {TT,CT}
    return CT
end

function Base.length(msum::MajoranaSumMulti{TT,CT}) where {TT<:Integer,CT}
    total_strings = 0
    for wi in eachindex(msum.MultiMajoranas)
        nstrings = length(msum.MultiMajoranas[wi])
        total_strings += nstrings
    end
    return total_strings
end

PropagationBase.storage(msum::MajoranaSumMulti) = msum.MultiMajoranas

function norm(msum::MajoranaSumMulti{TT,CT}, L=2) where {TT<:Integer,CT}
    if length(msum) == 0
        return 0.0
    end
    all_coeffs = zeros(CT, length(msum))
    idx = 1
    for wi in eachindex(msum.MultiMajoranas)
        dict = msum.MultiMajoranas[wi]
        length_of_sector = length(dict)
        all_coeffs[idx:idx+length_of_sector-1] .= collect(values(dict))
        idx += length_of_sector
    end
    return LinearAlgebra.norm(all_coeffs, L)
end

function Base.show(io::IO, msum::MajoranaSumMulti)
    max_display = 5
    print(io, "MajoranaSumMulti with $(length(msum)) terms\n")
    #=for (weight_key, dict) in msum.MultiMajoranas
        print(io, "Level $weight_key ($(length(dict)) terms):\n")
        for (i, (ms_int, coeff)) in enumerate(dict)
            print(io, "    $(coeff) * $(reverse(bitstring(ms_int)))\n")
            if i > max_display
                break 
            end
        end
    end=#
end

function show_stats(msum::MajoranaSumMulti{TT,CT}) where {TT<:Integer,CT}
    total_strings = 0
    for wi in eachindex(msum.MultiMajoranas)
        dict = msum.MultiMajoranas[wi]
        nstrings = length(dict)
        total_strings += nstrings
    end
    for wi in eachindex(msum.MultiMajoranas)
        dict = msum.MultiMajoranas[wi]
        nstrings = length(dict)
        println("Level $wi: $nstrings strings ($(round(100. * nstrings / total_strings; digits=2))%)")
    end
    println("Total strings: $total_strings")
end

function PropagationBase.nsites(msum::MajoranaSumMulti)
    return msum.nsites
end

function nfermions(msum::MajoranaSumMulti)
    if msum.is_spinful
        return 2 * msum.nsites
    else
        return msum.nsites
    end
end

function is_spinful(msum::MajoranaSumMulti)
    return msum.is_spinful
end
