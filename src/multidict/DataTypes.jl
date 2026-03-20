import Base: keys

struct MajoranaSumMulti{TT<:Integer,CT} <: AbstractMajoranaSum
    nsites::Int
    is_spinful::Bool
    MultiMajoranas::Dict{Int64,Dict{TT,CT}}
end

function MajoranaSumMulti(msum::MajoranaSum{TT,CT}, level_mapper::Function, n_levels::Int) where {TT<:Integer,CT}
    nsites = msum.nsites
    is_spinful = msum.is_spinful
    multimajs = Dict{Int64,Dict{TT,CT}}()

    for (ms_int, coeff) in msum.Majoranas
        ms_weight = get_weight(ms_int)
        ms_level = level_mapper(ms_int)
        ms_key = (ms_weight - 1) * n_levels + ms_level
        if !haskey(multimajs, ms_key)
            multimajs[ms_key] = Dict{TT,CT}()
        end
        multimajs[ms_key][ms_int] = coeff
    end
    return MajoranaSumMulti{TT,CT}(nsites, is_spinful, multimajs)
end

function Base.keys(msum::MajoranaSumMulti)
    return Base.keys(msum.MultiMajoranas)
end

# set ms of certain weight assuming dict weight is already present in msum
function set!(msum::MajoranaSumMulti{TT,CT}, weight_key, ms_int::TT, coeff::CT) where {TT<:Integer,CT}
    msum.MultiMajoranas[weight_key][ms_int] = coeff
end

function similar(msum::MajoranaSumMulti{TT,CT}, W::Int, nlevels::Int) where {TT<:Integer,CT}
    out_dict = Dict{Int64,Dict{TT,CT}}()
    for weight_key in (W - 2, W, W + 2)
        for j = 1:nlevels
            out_dict[(weight_key - 1) * nlevels + j] = Dict{TT,CT}()
        end
    end
    return MajoranaSumMulti(msum.nsites, msum.is_spinful, out_dict)
end

function coefftype(msum::MajoranaSumMulti{TT,CT}) where {TT,CT}
    return CT
end

function Base.length(msum::MajoranaSumMulti{TT,CT}) where {TT<:Integer,CT}
    total_strings = 0
    for (weight_key, dict) in msum.MultiMajoranas
        nstrings = length(dict)
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
    for (weight_key, dict) in msum.MultiMajoranas
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

function show_stats(msum::MajoranaSumMulti{TT,CT}, n_levels::Int) where {TT<:Integer,CT}
    total_strings = 0
    for (weight_key, dict) in msum.MultiMajoranas
        nstrings = length(dict)
        total_strings += nstrings
    end
    sorted_keys = sort(collect(keys(msum.MultiMajoranas)))
    for weight_key in sorted_keys
        weight_sector, level = _split_key(weight_key, n_levels)
        nstrings = length(msum.MultiMajoranas[weight_key])
        println("Level $weight_sector - $level: $nstrings strings ($(round(100. * nstrings / total_strings))%)")
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

function _split_key(key::Int64, n_levels::Int)
    weight = div(key - 1, n_levels) + 1
    level = mod(key - 1, n_levels) + 1
    return weight, level
end

function _get_weight_from_key(key::Int64, n_levels::Int)
    return _split_key(key, n_levels)[1]
end

function _get_level_from_key(key::Int64, n_levels::Int)
    return _split_key(key, n_levels)[2]
end


