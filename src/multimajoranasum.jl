###
##
# A Majorana sum split over work zones, each a `MajoranaSum` or `VectorMajoranaSum` that one thread owns for
# the whole propagation, as `MultiPauliSum` in PauliPropagation. PropagationBase assigns every Majorana string
# to a zone by a map that is linear in its bits, so a Majorana rotation sends the strings it creates from one
# zone to exactly one other, and every operation runs zone by zone through the `MultiSumStorage` trait: the
# gates need no code of their own.
##
###

"""
    MultiMajoranaSum(msum::AbstractMajoranaSum, n_zones::Integer=defaultnzones())

Split the Majorana sum `msum` over `n_zones` work zones, each a sum of the type of `msum` that one thread owns for the whole propagation.
Propagating it parallelizes over the zones, with one hand-off of new Majorana strings between the zones per gate.
`n_zones` must be a power of two and defaults to one zone per thread.
`VectorMajoranaSum` zones are faster than `MajoranaSum` ones.
Convert back with `MajoranaSum(msum)` or `VectorMajoranaSum(msum)`.
"""
struct MultiMajoranaSum{TS<:AbstractMajoranaSum,ZM<:ZoneMap} <: AbstractMajoranaSum
    nsites::Int
    zones::Vector{TS}
    zonemap::ZM

    # a zone is looked up by the index the map assigns without a bounds check
    function MultiMajoranaSum(nsites::Integer, zones::Vector{TS}, zonemap::ZM) where {TS<:AbstractMajoranaSum,ZM<:ZoneMap}
        if length(zones) != nzones(zonemap)
            throw(ArgumentError("got $(length(zones)) zones for a zone map of $(nzones(zonemap))."))
        end
        return new{TS,ZM}(nsites, zones, zonemap)
    end
end

MultiMajoranaSum(msum::AbstractMajoranaSum, n_zones::Integer=defaultnzones()) = add!(_emptyzones(msum, n_zones), msum)

# splitting an already split sum re-zones it, so the zones are seeded from a zone and not from it
MultiMajoranaSum(msum::MultiMajoranaSum, n_zones::Integer=defaultnzones()) = add!(_emptyzones(first(zones(msum)), n_zones), msum)

# `seed` fixes the type of the zones
_emptyzones(seed::AbstractMajoranaSum, n_zones::Integer) =
    MultiMajoranaSum(nsites(seed), [emptylike(seed) for _ in 1:n_zones], ZoneMap(termtype(seed), n_zones))

# The zones carry the storage trait, and the spinfulness, since PropagationBase builds a multi sum around new
# zones as `MultiMajoranaSum(nsites, zones, zonemap)`.
PropagationBase.storage(msum::MultiMajoranaSum) = msum.zones
PropagationBase.nsites(msum::MultiMajoranaSum) = msum.nsites
is_spinful(msum::MultiMajoranaSum) = is_spinful(first(msum.zones))

# `TS(prop_cache)` extracts the sum and converts it, so gathering the zones has to be reachable through `convert`
Base.convert(::Type{MajoranaSum}, msum::MultiMajoranaSum) = MajoranaSum(msum)
Base.convert(::Type{VectorMajoranaSum}, msum::MultiMajoranaSum) = VectorMajoranaSum(msum)

wrapcoefficients(msum::MultiMajoranaSum, ::Type{MProp}) where {MProp<:PathProperties} =
    MultiMajoranaSum(msum.nsites, [wrapcoefficients(zone, MProp) for zone in zones(msum)], zonemap(msum))
unwrapcoefficients(msum::MultiMajoranaSum) = MultiMajoranaSum(msum.nsites, map(unwrapcoefficients, zones(msum)), zonemap(msum))

function Base.show(io::IO, msum::MultiMajoranaSum)
    println(io, "MultiMajoranaSum of $(nameof(eltype(zones(msum)))) with $(length(msum)) terms over $(nzones(msum)) zones:")
    println(io, "  zone sizes: ", zonesizes(msum))

    # the zones are printed in turn, so the Majorana strings do not come out in any particular order
    for (i, (ms, coeff)) in enumerate(msum)
        if i > 20
            println(io, "  ...")
            break
        end
        println(io, "  ", coeff, " * ", reverse(bitstring(ms)))
    end
end

"""
    MultiMajoranaPropagationCache(msum::MultiMajoranaSum)

Propagation cache for a `MultiMajoranaSum`: the propagation cache of every zone, and one outbox per zone, itself a `MultiMajoranaSum`, in which a zone leaves the Majorana strings it creates for the zones that own them.
It carries no auxiliary sum of its own, since every zone brings one.
"""
struct MultiMajoranaPropagationCache{MS<:MultiMajoranaSum,ZC<:AbstractMajoranaPropagationCache} <: AbstractMajoranaPropagationCache
    msum::MS
    zonecaches::Vector{ZC}
    outboxes::Vector{MS}

    # one cache and one box per zone, looked up by the index the zone map assigns
    function MultiMajoranaPropagationCache(msum::MS, zonecaches::Vector{ZC}, outboxes::Vector{MS}) where {MS<:MultiMajoranaSum,ZC<:AbstractMajoranaPropagationCache}
        if !(length(zonecaches) == length(outboxes) == nzones(msum))
            throw(ArgumentError("got $(length(zonecaches)) zone caches and $(length(outboxes)) outboxes for $(nzones(msum)) zones."))
        end
        return new{MS,ZC}(msum, zonecaches, outboxes)
    end
end

function MultiMajoranaPropagationCache(msum::MultiMajoranaSum)
    zonecaches = map(PropagationCache, zones(msum))
    outboxes = [emptylike(msum) for _ in 1:nzones(msum)]
    return MultiMajoranaPropagationCache(msum, zonecaches, outboxes)
end

PropagationBase.PropagationCache(msum::MultiMajoranaSum) = MultiMajoranaPropagationCache(msum)

PropagationBase.mainsum(prop_cache::MultiMajoranaPropagationCache) = prop_cache.msum
PropagationBase.auxsum(prop_cache::MultiMajoranaPropagationCache) = throw(ArgumentError("MultiMajoranaPropagationCache does not have an auxiliary sum."))

function Base.show(io::IO, prop_cache::MultiMajoranaPropagationCache)
    println(io, "MultiMajoranaPropagationCache with $(length(prop_cache)) terms over $(nzones(prop_cache)) zones:")
    println(io, "  zone sizes: ", zonesizes(prop_cache))
end
