###
## Majorana Propagation Surrogate
##
# Mirrors the (experimental) Pauli propagation Surrogate of PauliPropagation.jl, but for
# Majorana propagation. Instead of propagating numerical coefficients, we propagate a symbolic
# graph of `cos`/`sin` nodes. The graph is built *once* by a single `propagate`, then evaluated
# cheaply for *many* parameter vectors via `evaluate!`. This is the analog of
# PauliPropagation's `NodePathProperties` surrogate and reuses as much of it as possible.
#
# We reuse from PauliPropagation:
#   - the leaf node type `EvalEndNode` and the abstract `CircuitNode`,
#   - the node accessors `parents`, `getnodeval`, `isevaluated`, `setcummulativevalue`,
#   - the recursive `reset!` over the node graph and `_traceevalorder(::EvalEndNode, ...)`.
#
# The only genuinely new piece is `MajoranaRotationNode`: it is the analog of
# `PauliRotationNode` but additionally carries a per-edge angle `scale`. This is required because
# a `FermionicRotation` decomposes into `MajoranaRotation`s that rotate by `theta * coeff * 2`
# (see `applymergetruncate!(::FermionicRotation, ...)` in `gates.jl`), and PauliPropagation's
# `PauliRotationNode`/`NodePathProperties` are a closed, scale-free union we cannot extend.
###

using PauliPropagation: CircuitNode, EvalEndNode, getnodeval, isevaluated, setcummulativevalue
import PauliPropagation: _traceevalorder, reset!, evaluate!

## ----------------------------------------------------------------------------------------------
## Node and coefficient types
## ----------------------------------------------------------------------------------------------

"""
    MajoranaRotationNode

Surrogate graph node for a `MajoranaRotation`. Analogous to PauliPropagation's `PauliRotationNode`
but carries a per-edge angle `scale`, so that a `FermionicRotation` (which rotates by
`theta * coeff * 2`) can be represented. For a plain `MajoranaRotation` the scale is `1`.
The trig factor of edge `i` is `cos`/`sin` of `scales[i] * thetas[param_idx]` (`trig_inds[i]` is
`1` for cos, `-1` for sin), times the integer `signs[i]`.
"""
mutable struct MajoranaRotationNode <: CircuitNode
    parents::Vector{Union{EvalEndNode,MajoranaRotationNode}}
    trig_inds::Vector{Int}        # 1 = cos, -1 = sin
    signs::Vector{Int}            # ±1 prefactor (from the ±i prefactor of `ms_mult`)
    scales::Vector{Float64}       # per-edge angle scale (coeff*2 for FermionicRotation, 1 for MajoranaRotation)
    param_idx::Int
    cummulative_value::Float64    # change to Real for automatic differentiation
    is_evaluated::Bool
end

const _MajParent = Union{EvalEndNode,MajoranaRotationNode}

Base.show(io::IO, node::MajoranaRotationNode) =
    print(io, "MajoranaRotationNode($(length(node.parents)) parent(s), param_idx=$(node.param_idx))")

"""
    MajoranaNodePathProperties(node)
    MajoranaNodePathProperties(node, nsins, ncos, freq)

Surrogate `PathProperties` type for Majorana propagation. Carries a `CircuitNode` instead of a
numerical coefficient. If `nsins`, `ncos` and `freq` are not given they are initialized to 0.
"""
struct MajoranaNodePathProperties <: PathProperties
    node::_MajParent
    nsins::Int
    ncos::Int
    freq::Int
end

MajoranaNodePathProperties(node::_MajParent) = MajoranaNodePathProperties(node, 0, 0, 0)

PropagationBase.numcoefftype(::Type{MajoranaNodePathProperties}) = Float64

"""
    tonumber(path::MajoranaNodePathProperties)

Get the cumulative coefficient of a `MajoranaNodePathProperties`. Assumes the surrogate has
already been evaluated with `evaluate!`.
"""
PropagationBase.tonumber(path::MajoranaNodePathProperties) = path.node.cummulative_value

Base.show(io::IO, pth::MajoranaNodePathProperties) =
    print(io, "MajoranaNodePathProperties($(typeof(pth.node)), nsins=$(pth.nsins), ncos=$(pth.ncos), freq=$(pth.freq))")

## ----------------------------------------------------------------------------------------------
## Wrapping coefficients
## ----------------------------------------------------------------------------------------------

"""
    wrapcoefficients(msum::MajoranaSum, ::Type{MajoranaNodePathProperties})

Wrap the coefficients of a `MajoranaSum` into `MajoranaNodePathProperties`, turning it into the
seed of a surrogate. Coefficients are assumed real (stored as `Float64`), mirroring the
PauliPropagation surrogate.
"""
function wrapcoefficients(msum::MajoranaSum, ::Type{MajoranaNodePathProperties})
    if length(msum) == 0
        throw("The MajoranaSum is empty.")
    end
    # `EvalEndNode.pstr::Int` is informational only (used for `show`); store the Majorana string
    # integer best-effort. It is never read during evaluation/merge/truncation.
    terms = Dict(ms => MajoranaNodePathProperties(EvalEndNode(ms % Int, float(coeff), 0.0, false))
                 for (ms, coeff) in msum)
    return MajoranaSum(nsites(msum), is_spinful(msum), terms)
end

## ----------------------------------------------------------------------------------------------
## Node builders and cos/sin application (analog of Surrogate/propagate.jl `_applycos`/`_applysin`)
## ----------------------------------------------------------------------------------------------

# Build a fresh rotation node with a single edge.
function _buildrotationnode(node::_MajParent, trig_ind::Int, sign::Int, scale::Real, param_idx::Int)
    return MajoranaRotationNode(_MajParent[node], Int[trig_ind], Int[sign], Float64[scale], param_idx, 0.0, false)
end

# `theta` carries the dummy parameter index; `scale` is the per-edge angle scale.
function _applycos(path::MajoranaNodePathProperties, param_idx::Integer, sign::Integer=1, scale::Real=1.0)
    node = _buildrotationnode(path.node, 1, Int(sign), scale, Int(param_idx))
    return MajoranaNodePathProperties(node, path.nsins, path.ncos + 1, path.freq + 1)
end

function _applysin(path::MajoranaNodePathProperties, param_idx::Integer, sign::Integer=1, scale::Real=1.0)
    node = _buildrotationnode(path.node, -1, Int(sign), scale, Int(param_idx))
    return MajoranaNodePathProperties(node, path.nsins + 1, path.ncos, path.freq + 1)
end

## ----------------------------------------------------------------------------------------------
## Merging (analog of Surrogate/propagate.jl `mergefunc`/`_mergenodes!`)
## ----------------------------------------------------------------------------------------------

function _mergenodes!(node1::MajoranaRotationNode, node2::MajoranaRotationNode)
    append!(node1.parents, node2.parents)
    append!(node1.trig_inds, node2.trig_inds)
    append!(node1.signs, node2.signs)
    append!(node1.scales, node2.scales)
    return node1
end

function PropagationBase.mergefunc(pth1::MajoranaNodePathProperties, pth2::MajoranaNodePathProperties)
    return MajoranaNodePathProperties(
        _mergenodes!(pth1.node, pth2.node),
        min(pth1.nsins, pth2.nsins),
        min(pth1.ncos, pth2.ncos),
        min(pth1.freq, pth2.freq),
    )
end

## ----------------------------------------------------------------------------------------------
## Gate application
## ----------------------------------------------------------------------------------------------

const MajSurrogateCache = MajoranaPropagationCache{<:MajoranaSum{<:Integer,MajoranaNodePathProperties}}

"""
Surrogate action of a `MajoranaRotation`. Instead of multiplying by `cos`/`sin` values, it builds
`cos` and `sin` graph nodes that record the parameter index `theta` and the angle `scale`.
Mirrors the numeric `applytoall!(::MajoranaRotation, ...)` in `gates.jl`.
"""
function PropagationBase.applytoall!(gate::MajoranaRotation, prop_cache::MajSurrogateCache, theta; scale=1.0, kwargs...)
    msum = mainsum(prop_cache)
    aux_msum = auxsum(prop_cache)

    param_idx = theta  # the dummy thetas transport the parameter indices
    gate_int = gate.ms_int

    for (ms_int, coeff) in msum
        if commutes(gate_int, ms_int)
            continue
        end

        # exp(i theta gate/2) ms exp(-i theta gate/2) = cos(theta) ms + i sin(theta) gate*ms
        prefactor, new_ms = ms_mult(gate_int, ms_int, nfermions(msum))

        # cos branch keeps the original Majorana string
        set!(msum, ms_int, _applycos(coeff, param_idx, 1, scale))
        # sin branch creates the new (non-overlapping) Majorana string; sign = -imag(±i) = ±1
        set!(aux_msum, new_ms, _applysin(coeff, param_idx, Int(-imag(prefactor)), scale))
    end

    return
end

"""
Surrogate action of a `FermionicRotation`. Mirrors the numeric
`applymergetruncate!(::FermionicRotation, ...)` in `gates.jl`, but passes the per-sub-rotation
angle as a `scale` (= `coeff * 2`) while keeping `theta` as the dummy parameter index.
"""
function PropagationBase.applymergetruncate!(gate::FermionicRotation, prop_cache::MajSurrogateCache, theta; truncate_each_mr=nothing, kwargs...)
    ms_rotations, coeffs, truncate_after_each_majrot = getmajoranarotations(gate, nsites(prop_cache))
    if !isnothing(truncate_each_mr)
        truncate_after_each_majrot = truncate_each_mr
    end

    for (gate_ms, coeff) in zip(ms_rotations, coeffs)
        # `MajoranaRotation` implements exp(-i theta/2 mstring), hence the factor of 2
        applytoall!(gate_ms, prop_cache, theta; scale=coeff * 2.0, kwargs...)
        merge!(prop_cache; kwargs...)
        if truncate_after_each_majrot
            truncate!(prop_cache; kwargs...)
        end
    end
    if !truncate_after_each_majrot
        truncate!(prop_cache; kwargs...)
    end

    return prop_cache
end

## ----------------------------------------------------------------------------------------------
## Propagation entry points (analog of Surrogate/propagate.jl)
## ----------------------------------------------------------------------------------------------

function _checkmajoranasurrogationconditions(circ)
    if !all(isa(gate, MajoranaRotation) || isa(gate, FermionicRotation) for gate in circ)
        throw(ArgumentError("The Majorana surrogate currently only accepts `MajoranaRotation`s and `FermionicRotation`s."))
    end
    return
end

"""
    propagate(circ, msum::MajoranaSum{<:Integer,MajoranaNodePathProperties}; kwargs...)

Construct a Majorana propagation surrogate of `msum` propagated through `circ` in the Heisenberg
picture. The circuit must only contain `MajoranaRotation`s and `FermionicRotation`s. Truncations
based on a numerical coefficient value cannot be used. Everything else matches `propagate!()` for
the non-surrogate code. Evaluate the resulting surrogate with `evaluate!`.
"""
function PropagationBase.propagate(circ, msum::MajoranaSum{TT,MajoranaNodePathProperties}; kwargs...) where {TT<:Integer}
    _checkmajoranasurrogationconditions(circ)
    return propagate!(circ, deepcopy(msum); kwargs...)
end

"""
    propagate!(circ, msum::MajoranaSum{<:Integer,MajoranaNodePathProperties}; kwargs...)

In-place construction of a Majorana propagation surrogate (see [`propagate`](@ref)).
"""
function PropagationBase.propagate!(circ, msum::MajoranaSum{TT,MajoranaNodePathProperties}; kwargs...) where {TT<:Integer}
    _checkmajoranasurrogationconditions(circ)
    # dummy parameters transport the parameter indices through the standard propagation flow
    dummy_thetas = collect(Int, 1:countparameters(circ))
    return propagate!(circ, msum, dummy_thetas; kwargs...)
end

## ----------------------------------------------------------------------------------------------
## Evaluation (reuses PauliPropagation's `reset!` and `_traceevalorder(::EvalEndNode, ...)`)
## ----------------------------------------------------------------------------------------------

# Trace the cumulative value of a `MajoranaRotationNode` by recursively evaluating its parents.
# Reuses PauliPropagation's node accessors and the `EvalEndNode` leaf method.
function _traceevalorder(node::MajoranaRotationNode, thetas; eval_list=nothing)
    if isevaluated(node)
        return node.cummulative_value
    end

    val = 0.0
    for ii in eachindex(node.parents)
        angle = node.scales[ii] * thetas[node.param_idx]
        trig = node.trig_inds[ii] == 1 ? cos(angle) :
               node.trig_inds[ii] == -1 ? sin(angle) : 1.0
        this_val = node.signs[ii] * trig

        parent = node.parents[ii]
        other_val = isevaluated(parent) ? getnodeval(parent) : _traceevalorder(parent, thetas; eval_list)
        val += this_val * other_val
    end

    setcummulativevalue(node, val)
    node.is_evaluated = true
    if !isnothing(eval_list)
        push!(eval_list, node)
    end
    return node.cummulative_value
end

"""
    evaluate!(msum::MajoranaSum{<:Integer,MajoranaNodePathProperties}, thetas; reset=true)

Evaluate a Majorana surrogate at the parameter vector `thetas` (given in Schrödinger / circuit
order). After this call the cumulative values are stored on the nodes and can be read with
`tonumber`, e.g. via `overlapwithfock(msum, fock_state)`. If `reset=false`, the `is_evaluated`
flags are not cleared first; only do this if they have been reset manually.
"""
function evaluate!(msum::MajoranaSum{TT,MajoranaNodePathProperties}, thetas; reset=true) where {TT<:Integer}
    paths = collect(coefficients(msum))

    if reset
        for pth in paths
            reset!(pth.node)
        end
    end

    for pth in paths
        _traceevalorder(pth.node, thetas)
    end

    return msum
end

"""
    reset!(msum::MajoranaSum{<:Integer,MajoranaNodePathProperties})

Reset the nodes of a Majorana surrogate. Needs to be done between evaluations with different
parameters (handled automatically by `evaluate!` unless `reset=false`).
"""
function reset!(msum::MajoranaSum{TT,MajoranaNodePathProperties}) where {TT<:Integer}
    for pth in collect(coefficients(msum))
        reset!(pth.node)
    end
    return
end
