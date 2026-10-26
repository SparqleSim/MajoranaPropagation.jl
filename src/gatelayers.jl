###
##
# Layers of commuting Majorana and fermionic rotations, applied in one pass over the sum, class by class, as
# PropagationBase applies a `GateLayer` of commuting gates (its Layers/classes.jl). A rotation acts on the Majorana
# strings that share an odd number of Majoranas with its even-weight generator, and flips the bits of the generator on
# them, so the generators are both the masks of the rotations and what PropagationBase's `PrecomputedLookup` finds the
# acting rotations of a string from. A fermionic rotation enters the pass as its Majorana rotations and truncates as it
# does on its own: after each of them, or, for the hoppings, only after its last. Which gates of a circuit form a layer
# is decided by `tolayers` of PauliPropagation, from `qinds` and `commutes` below.
##
###

# The Majoranas a gate acts on, counted from one, which `GateLayer` and `tolayers` compare: gates on different Majoranas
# commute, since their generators have even weight.
PropagationBase.qinds(gate::Union{MajoranaRotation,ImaginaryMajoranaRotation}) = _majoranaindices(gate.ms_int)
PropagationBase.qinds(gate::Union{FermionicRotation,ImaginaryFermionicRotation}) = reduce(union, _generatorindices(gate); init=Int[])

# the Majoranas of every generator of a gate
_generatorindices(gate::MajoranaRotation) = [_majoranaindices(gate.ms_int)]
_generatorindices(gate::Union{FermionicRotation,ImaginaryFermionicRotation}) =
    [_majoranaindices(ms) for (ms, _) in MajoranaSum(maximum(gate.sites), gate.symbol, gate.sites)]

_majoranaindices(bits) = PropagationBase._bitpositions(bits) .+ 1

# Two rotations commute if every generator of one shares an even number of Majoranas with every generator of the other,
# which `tolayers` asks of gates that share Majoranas.
PauliPropagation.commutes(gate_a::Union{MajoranaRotation,FermionicRotation}, gate_b::Union{MajoranaRotation,FermionicRotation}) =
    all(iseven(length(intersect(a, b))) for a in _generatorindices(gate_a), b in _generatorindices(gate_b))

"""
    applymergetruncate!(layer::GateLayer, prop_cache::AbstractMajoranaPropagationCache, thetas; truncate_each_mr=nothing, thread=true, kwargs...)

Apply a layer of commuting `MajoranaRotation`s and `FermionicRotation`s in one pass over the Majorana sum, class by class, instead of one rotation after the other.
A class is a set of Majorana strings that the rotations of the layer map into each other; it is gathered into a small sum of its own, to which the rotations acting on it are applied in the order of the layer.
Every fermionic rotation truncates as it does on its own (see `truncate_each_mr`), so the result is that of the gates one after the other, except that strings no rotation acts on are truncated too.
Layers are built with `GateLayer` or gathered from a circuit by `tolayers`, which does not reorder gates: a layer pays off when it holds many rotations, so list commuting gates next to each other, such as the hoppings of disjoint bonds.
Layers of other gates, layers whose Majorana rotations do not all commute, and coefficients that are not numbers, such as `MajoranaFrequencyTracker`, apply their gates one after the other.
"""
function PropagationBase.applymergetruncate!(layer::GateLayer, prop_cache::AbstractMajoranaPropagationCache, thetas; truncate_each_mr=nothing, thread::Bool=true, kwargs...)
    plan = _classplan(layer, prop_cache, thetas, truncate_each_mr)
    if isnothing(plan)
        PropagationBase._applygatesonebyone!(layer, prop_cache, thetas; truncate_each_mr, thread, kwargs...)
        return prop_cache
    end
    truncfunc, _ = _truncfuncs(prop_cache; thread, kwargs...)
    _applylayerpass!(prop_cache, plan..., truncfunc; thread)
    return prop_cache
end

# The plan of a layer applied class by class, as `(plan, gate_of, truncates_each)`, or `nothing` for a layer applied gate
# by gate: every gate must be a Majorana or fermionic rotation, the sum one that PropagationBase groups by class, with
# numbers as coefficients, and the Majorana rotations must commute. Rotation `k` of the plan comes from gate `gate_of[k]`
# of the layer, which truncates after each of its rotations if `truncates_each[gate_of[k]]`, and after its last otherwise.
function _classplan(layer::GateLayer, prop_cache, thetas, truncate_each_mr)
    if !(all(gate -> gate isa Union{MajoranaRotation,FermionicRotation}, layer.gates) && PropagationBase._canapplylayer(prop_cache))
        return nothing
    end

    rotations = MajoranaRotation[]
    angles = []
    gate_of = Int[]
    truncates_each = Bool[]
    function pushrotations!(gate::MajoranaRotation, theta)
        push!(truncates_each, true)
        push!(rotations, gate)
        push!(angles, theta)
        push!(gate_of, length(truncates_each))
        return
    end
    function pushrotations!(gate::FermionicRotation, theta)
        gate_rotations, coeffs, truncate_after_each_majrot = getmajoranarotations(gate, nsites(prop_cache))
        push!(truncates_each, something(truncate_each_mr, truncate_after_each_majrot))
        for (rotation, coeff) in zip(gate_rotations, coeffs)
            push!(rotations, rotation)
            # the angle that `applymergetruncate!` gives the rotation of a fermionic gate
            push!(angles, theta * coeff * 2.0)
            push!(gate_of, length(truncates_each))
        end
        return
    end
    PropagationBase._foreachgate(pushrotations!, layer, thetas)

    masks = [_branchmask(rotation, prop_cache) for rotation in rotations]
    # the classes are those of commuting rotations
    for i in eachindex(masks), j in i+1:lastindex(masks)
        if !_commutes_evengate(masks[i], masks[j])
            return nothing
        end
    end
    rules = [_branchrule(rotation, prop_cache, angle) for (rotation, angle) in zip(rotations, angles)]
    plan = PropagationBase.LayerPlan(masks, rules, PropagationBase.PrecomputedLookup(masks))
    return plan, gate_of, truncates_each
end

# One pass of the layer over the sum, as PropagationBase's `_applypass!` makes it, with every class applied to by
# `_applytoclass!` in the scratch of its task.
function _applylayerpass!(prop_cache, plan, gate_of, truncates_each, truncfunc::F; thread::Bool) where {F}
    if length(prop_cache) == 0
        return prop_cache
    end
    TT = termtype(prop_cache)
    CT = coefftype(prop_cache)
    workspace = PropagationBase._takeworkspace(TT, CT, PropagationBase.ClassScratch{TT,CT})
    try
        sources = PropagationBase._passsources(StorageType(prop_cache), prop_cache, workspace, thread)
        representativeof(term) = PropagationBase._classrepresentative(plan, term)
        applytoclass!(sink, scratch, class_terms, class_coeffs) =
            _applytoclass!(sink, scratch.cache, plan, gate_of, truncates_each, truncfunc, class_terms, class_coeffs)
        PropagationBase._applytoclasses!(representativeof, applytoclass!, prop_cache, workspace, sources, thread)
    finally
        PropagationBase._putbackworkspace!(workspace)
    end
    return prop_cache
end

# The rotations of the plan that act on one class, `class_terms` with `class_coeffs`, applied as the gates of the layer
# apply them one by one: the class is loaded into the emptied sum of `cache`, every acting rotation `k` branches it by
# its rule, and its gate truncates the class after it if `truncates_each[gate_of[k]]`, or else once its rotations are
# done, before the next gate that acts on the class. What is left is written to `sink`. As in PropagationBase, a class
# that no rotation acts on is truncated as well.
function _applytoclass!(sink, cache, plan, gate_of, truncates_each, truncfunc::F, class_terms, class_coeffs) where {F}
    positions = PropagationBase._branchinggates(plan.lookup, first(class_terms))
    class_sum = mainsum(cache)
    empty!(cache)
    # the strings of a class differ in their pivot bits
    PropagationBase._distinguishby!(class_sum, plan.pivots)
    # room for twice the strings, which loading the class and branching it by its first rotation need at most
    resize!(cache, 2 * length(class_terms))
    for index in eachindex(class_terms, class_coeffs)
        add!(class_sum, class_terms[index], class_coeffs[index])
    end

    # the gate whose rotations acted on the class since it was last truncated, or 0
    pending = 0
    for word in eachindex(positions)
        bits = positions[word]
        while bits != 0
            k = 64 * (word - 1) + trailing_zeros(bits) + 1
            bits &= bits - one(UInt64)
            gate = gate_of[k]
            if pending != 0 && pending != gate
                mergeandtruncate!(truncfunc, cache; thread=false)
            end
            if truncates_each[gate]
                xorbranchmergeandtruncate!(plan.rules[k], truncfunc, cache, plan.gate_masks[k]; thread=false)
                pending = 0
            else
                xorbranchmergeandtruncate!(plan.rules[k], nothing, cache, plan.gate_masks[k]; thread=false)
                pending = gate
            end
        end
    end
    if pending != 0 || all(iszero, positions)
        mergeandtruncate!(truncfunc, cache; thread=false)
    end

    for (term, coeff) in mainsum(cache)
        PropagationBase._emit!(sink, term, coeff)
    end
    return sink
end
