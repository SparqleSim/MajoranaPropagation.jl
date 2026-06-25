"""
    MajoranaTransferMap{TT<:Integer,CT}

A contiguous lookup table ("transfer map") for the conjugate action of a fermionic gate
on Majorana strings, analogous to `PauliPropagation.TransferMap`.

The map is indexed by the *restricted* Majorana string on the gate's modes, compressed to
canonical bit positions `0 .. nmodes-1` (so `tmap[0]` is the image of the identity column).

Each column is a list of `(cumulative_string, coeff)` entries. Unlike the Pauli case, the
sign of multiplying two Majorana strings depends on the *full* string and not only on the
restricted part. Therefore an entry stores the *cumulative gate string* `cumulative_string`
(the Majorana string that, multiplied onto the observable, produces the output string) rather
than the output string and a fixed coefficient. At application time the output string and the
full-string-dependent sign are obtained from `ms_mult(cumulative_string, full_observable)`,
while `coeff` carries the remaining (full-string-independent) prefactor.

Storage is a "structure of arrays": `entries` holds all `(cumulative_string, coeff)` tuples
contiguously, and `offsets` indexes where each column starts (`offsets[i] .. offsets[i+1]-1`
are the entries of the 0-based column `i-1`).
"""
struct MajoranaTransferMap{TT<:Integer,CT}
    entries::Vector{Tuple{TT,CT}}
    offsets::Vector{Int}
end

"""
    MajoranaTransferMap(columns::Vector{Vector{Tuple{TT,CT}}})

Build a `MajoranaTransferMap` from a dense vector of columns, where `columns[i]` holds the
entries of the 0-based column `i-1`.
"""
function MajoranaTransferMap(columns::Vector{Vector{Tuple{TT,CT}}}) where {TT<:Integer,CT}
    total_entries = sum(length, columns; init=0)
    entries = Vector{Tuple{TT,CT}}(undef, total_entries)
    offsets = Vector{Int}(undef, length(columns) + 1)

    next_entry = 1
    for (column_index, column) in enumerate(columns)
        offsets[column_index] = next_entry
        for item in column
            entries[next_entry] = item
            next_entry += 1
        end
    end
    offsets[end] = next_entry

    return MajoranaTransferMap{TT,CT}(entries, offsets)
end

ncolumns(tmap::MajoranaTransferMap) = length(tmap.offsets) - 1
Base.length(tmap::MajoranaTransferMap) = length(tmap.entries)

"""
    getindex(tmap::MajoranaTransferMap, column_index::Integer)

Return a view of the entries of the 0-based `column_index`.
"""
function Base.getindex(tmap::MajoranaTransferMap, column_index::Integer)
    index = Int(column_index)
    if index < 0 || index >= ncolumns(tmap)
        throw(BoundsError(tmap, column_index))
    end
    start = tmap.offsets[index+1]
    stop = tmap.offsets[index+2] - 1
    return @view tmap.entries[start:stop]
end

function Base.:(==)(left::MajoranaTransferMap, right::MajoranaTransferMap)
    return left.entries == right.entries && left.offsets == right.offsets
end

function Base.show(io::IO, tmap::MajoranaTransferMap)
    print(io, "MajoranaTransferMap($(ncolumns(tmap)) columns, $(length(tmap)) entries)")
end

"""
    SurrogateCoeff(path::MajoranaNodePathProperties, inv_mu::ComplexF64)

Symbolic coefficient stored in an angle-free (surrogate-backed) transfer map. `path` is a
surrogate node graph encoding the angle dependence of the propagated coefficient `a`, while
`inv_mu = 1/μ` is the angle-independent frame prefactor. Evaluating the table at an angle gives the
numeric coefficient `ctilde = a * inv_mu` (see [`evaluate`](@ref)).
"""
struct SurrogateCoeff
    path::MajoranaNodePathProperties
    inv_mu::ComplexF64
end
