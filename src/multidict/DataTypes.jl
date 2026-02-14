import Base: keys

struct MajoranaSumMulti{TT<:Integer,CT} <: AbstractMajoranaSum
    nsites::Int
    is_spinful::Bool
    MultiMajoranas::Dict{Int,Dict{TT,CT}}
end

function MajoranaSumMulti(msum::MajoranaSum{TT,CT}) where {TT<:Integer,CT}
    nsites = msum.nsites
    is_spinful = msum.is_spinful
    multimajs = Dict{Int,Dict{TT,CT}}()

    for (ms_int, coeff) in msum.Majoranas
        ms_weight = get_weight(ms_int)
        if !haskey(multimajs, ms_weight)
            multimajs[ms_weight] = Dict{TT,CT}()
        end
        multimajs[ms_weight][ms_int] = coeff
    end
    return MajoranaSumMulti{TT,CT}(nsites, is_spinful, multimajs)
end

function Base.keys(msum::MajoranaSumMulti)
    return Base.keys(msum.MultiMajoranas)
end

# set ms of certain weight assuming dict weight is already present in msum
function set!(msum::MajoranaSumMulti{TT,CT}, weight_key::Int, ms_int::TT, coeff::CT) where {TT<:Integer,CT}
    msum.MultiMajoranas[weight_key][ms_int] = coeff
end

function similar(msum::MajoranaSumMulti{TT,CT}, W) where {TT<:Integer,CT}
    return MajoranaSumMulti(msum.nsites, msum.is_spinful,
        Dict(W - 2 => Dict{TT,CT}(),
            W => Dict{TT,CT}(),
            W + 2 => Dict{TT,CT}()))
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

function show_stats(msum::MajoranaSumMulti{TT,CT}) where {TT<:Integer,CT}
    total_strings = 0
    for (weight_key, dict) in msum.MultiMajoranas
        nstrings = length(dict)
        total_strings += nstrings
    end
    sorted_keys = sort(collect(keys(msum.MultiMajoranas)))
    for weight_key in sorted_keys
        nstrings = length(msum.MultiMajoranas[weight_key])
        println("Weight $weight_key: $nstrings strings ($(round(100. * nstrings / total_strings))%)")
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
