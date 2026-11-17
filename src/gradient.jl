###
##
# Gradients of an expectation value by the angles of a circuit, from one forward and one backward sweep, as
# `rewindgradient` of PauliPropagation computes them (its gradient.jl). The backward sweep, its state and the
# coefficient that carries the dual beside the coefficient of the operator are PauliPropagation's; what is here is
# what depends on the Majorana algebra: the dual of a Majorana string, the gradient of a Majorana rotation, and how
# Majorana and fermionic rotations are undone.
##
###

using PauliPropagation: _DualCoeff, _BackwardSweepState, _KeepAbsDual, _zerogradient, _layergradient, _undogatesonebyone!, GradientScratch

"""
    rewindgradient(circuit, msum::AbstractMajoranaSum, thetas, overlapfunc; thread=true, keep_rel_dual=0.03, kwargs...)

Compute `overlapfunc(propagate(circuit, msum, thetas; kwargs...))` together with its gradient with respect to `thetas`, in one forward and one backward sweep costing ``O(\\mathrm{length}(circuit))`` gate applications.
Returns `(expec, grad)`, where `grad` has the shape of `thetas`, with a vector of derivatives for every `GateLayer`.
The circuit may hold `MajoranaRotation`s, `FermionicRotation`s and `GateLayer`s of them, and the coefficients of `msum` must be numbers.
`overlapfunc` must be linear in the coefficients of the Majorana sum it is handed, e.g. `msum -> overlapwithfock(msum, fock_state)` or `msum -> scalarproduct(msum, target)`.
`kwargs` are the truncations of `propagate`, which apply to both sweeps, so the gradient is that of the truncated propagation only approximately.
The backward sweep keeps a Majorana string that `min_abs_coeff` or `min_rel_coeff` would drop as long as the derivative of the expectation by its coefficient is, in absolute value, at least `keep_rel_dual` times the largest such derivative over the strings of the propagated operator; `keep_rel_dual=Inf` truncates by the coefficients of the operator alone.

This design is adapted from the publication "Backpropagating Pauli Propagation" by Lin et al. (arXiv:2607.15184).
"""
function PauliPropagation.rewindgradient(circuit, msum::AbstractMajoranaSum, thetas, overlapfunc; kwargs...)
    return rewindgradient!(circuit, deepcopy(msum), thetas, overlapfunc; kwargs...)
end

"""
    rewindgradient!(circuit, msum::AbstractMajoranaSum, thetas, overlapfunc; kwargs...)
    rewindgradient!(circuit, prop_cache::AbstractMajoranaPropagationCache, thetas, overlapfunc; kwargs...)

In-place version of `rewindgradient`, which propagates `msum` or `prop_cache` itself and leaves the propagated operator in `prop_cache`.
"""
function PauliPropagation.rewindgradient!(circuit, msum::AbstractMajoranaSum, thetas, overlapfunc; kwargs...)
    return rewindgradient!(circuit, PropagationCache(msum), thetas, overlapfunc; kwargs...)
end

function PauliPropagation.rewindgradient!(circuit, forward_cache::AbstractMajoranaPropagationCache, thetas, overlapfunc;
    thread::Bool=true, keep_rel_dual::Real=0.03, customtruncfunc=nothing, kwargs...)

    circuit, thetas = PropagationBase._promotecircandparams(circuit, thetas)
    if !all(_isrewindable, circuit)
        throw(ArgumentError("rewindgradient undoes MajoranaRotations, FermionicRotations and GateLayers of them only."))
    end
    if !(coefftype(forward_cache) <: Number)
        throw(ArgumentError("rewindgradient needs numbers as coefficients, got $(coefftype(forward_cache))."))
    end

    # forward sweep: ordinary Heisenberg propagation, exactly as in `propagate`
    propagate!(circuit, forward_cache, thetas; thread, customtruncfunc, kwargs...)
    expec = overlapfunc(activesum(forward_cache))

    # the dual sum starts from the final operator, with overlapfunc applied to each of its Majorana strings individually,
    # and is carried beside the coefficients of the operator through the backward sweep
    cache = PauliPropagation.wrapdual(forward_cache, overlapfunc; thread)
    state = _BackwardSweepState(cache, _zerogradient(thetas), 0)

    # without a nonzero dual coefficient the gradient is zero
    max_abs_dual = mapreducecoeffs(coeff -> abs(coeff.dual), max, cache; thread)
    if iszero(max_abs_dual)
        return expec, state.grad
    end

    # backward sweep: the gates are undone in the order of the circuit, the reverse of the Heisenberg sweep, each by
    # its own angle, which the gate negates as it is undone
    keep = _KeepAbsDual(keep_rel_dual * max_abs_dual, customtruncfunc)
    PropagationBase._propagate!(PauliPropagation._undostep!, circuit, state, thetas; thread, customtruncfunc=keep, kwargs...)

    return expec, state.grad
end

# the gates that the backward sweep undoes and differentiates by
_isrewindable(gate) = false
_isrewindable(gate::Union{MajoranaRotation,FermionicRotation}) = true
_isrewindable(layer::GateLayer) = all(_isrewindable, layer.gates)

# The truncation of coefficients in the backward sweep: a string below the cutoff is dropped only if its dual
# coefficient is below `keep_abs_dual` as well, and the custom truncation of the user reads the coefficient of the operator.
_belowcutoff(keep::_KeepAbsDual, coeff, cutoff) =
    cutoff > 0 && PauliPropagation.truncatemincoeff(coeff, cutoff) && PauliPropagation.truncatemincoeff(coeff.dual, keep.keep_abs_dual)
_customtruncates(keep::_KeepAbsDual, mstr, coeff) = _customtruncates(keep.customtruncfunc, mstr, coeff.coeff)


# The operator of `prop_cache` with its dual: every Majorana string keeps its coefficient and gets the overlap it has
# on its own as its coefficient in the dual sum.
function PauliPropagation.wrapdual(prop_cache::AbstractMajoranaPropagationCache, overlapfunc; thread::Bool=true)
    op_sum = activesum(prop_cache)
    n_sites = nsites(op_sum)
    spinful = is_spinful(op_sum)
    # ponytail: one small sum per string, since an overlap is a function of a sum; give the overlaps a method for a
    # single string if building the dual shows up in a profile
    singletonoverlap(term) = overlapfunc(VectorMajoranaSum(n_sites, spinful, [term], [1.0], 1))

    # the type of an overlap, from that of the identity
    dual_type = float(typeof(singletonoverlap(zero(termtype(prop_cache)))))
    coeff_type = _DualCoeff{coefftype(prop_cache),dual_type}
    cache = PropagationCache(PauliPropagation.convertcoefftype(coeff_type, op_sum))
    withoverlap(term, coeff) = coeff_type(coeff.coeff, singletonoverlap(term))
    mapcoeffsbypair!(withoverlap, cache; thread)
    return cache
end


# Records the gradient of a Majorana rotation and undoes it.
function PauliPropagation._undostep!(gate::MajoranaRotation, state::_BackwardSweepState, theta; thread::Bool=true, kwargs...)
    state.param_idx += 1
    state.grad[state.param_idx] = _generatorcommutatordot(gate, state.cache; thread)

    # the operator sum truncates normally, and the dual sum loses every Majorana string that it loses
    applymergetruncate!(gate, state.cache, -theta; thread, kwargs...)
    return state
end

# Records the gradient of a fermionic rotation and undoes it: its Majorana rotations are undone in reverse, each by
# its negated angle, and truncate as `applymergetruncate!` has them truncate. The angle of the gate enters every one
# of them, with twice its coefficient.
function PauliPropagation._undostep!(gate::FermionicRotation, state::_BackwardSweepState, theta; truncate_each_mr=nothing,
    thread::Bool=true, kwargs...)

    cache = state.cache
    ms_rotations, coeffs, truncate_after_each_majrot = getmajoranarotations(gate, nsites(cache))
    truncate_after_each_majrot = something(truncate_each_mr, truncate_after_each_majrot)
    truncfunc, termdrop = _truncfuncs(cache; thread, kwargs...)

    grad = 0.0
    for k in reverse(eachindex(ms_rotations))
        is_last = k == firstindex(ms_rotations)
        truncates = truncate_after_each_majrot || is_last
        grad += 2.0 * coeffs[k] * _generatorcommutatordot(ms_rotations[k], cache; thread)
        _xorbranchstep!(ms_rotations[k], cache, -theta * coeffs[k] * 2.0, truncates ? truncfunc : nothing, truncates ? termdrop : nothing;
            force_truncate=is_last && !truncate_after_each_majrot, thread)
    end

    state.param_idx += 1
    state.grad[state.param_idx] = grad
    return state
end

# Records the gradient of a layer, the derivatives by the angles of its gates, and undoes it as `applymergetruncate!`
# applies it: class by class, with the negated angles, or else its gates one by one, in the order of the layer, since
# the Heisenberg sweep applied them in reverse.
function PauliPropagation._undostep!(layer::GateLayer, state::_BackwardSweepState{<:AbstractMajoranaPropagationCache}, thetas;
    truncate_each_mr=nothing, min_rel_coeff=nothing, thread::Bool=true, kwargs...)

    cache = state.cache
    plan = isnothing(min_rel_coeff) ? _classplan(layer, cache, -thetas, truncate_each_mr) : nothing
    gate_grads = if isnothing(plan)
        _undogatesonebyone!(layer, cache, thetas; truncate_each_mr, min_rel_coeff, thread, kwargs...)
    else
        truncfunc, _ = _truncfuncs(cache; thread, kwargs...)
        _undolayerpass!(layer, cache, thetas, plan..., truncfunc; thread)
    end

    state.param_idx += 1
    state.grad[state.param_idx] = _layergradient(thetas, gate_grads)
    return state
end

# One pass of the layer over the sum as `_applylayerpass!` makes it, in which every task applies the rotations through
# rules of its own, PauliPropagation's `_GradientRule`s, which add what every two strings a rotation mixes contribute
# to the gradient of that rotation. The rotations of a class plan commute, so each of them may record its gradient
# whenever the pass gets to it. Returns the gradients of the gates of the layer.
function _undolayerpass!(layer, cache, thetas, plan, gate_of, truncates_each, truncfunc::F; thread::Bool) where {F}
    gate_grads = zeros(length(layer.gates))
    if length(cache) == 0
        return gate_grads
    end

    n_rotations = length(plan.rules)
    TT = termtype(cache)
    CT = coefftype(cache)
    workspace = PropagationBase._takeworkspace(TT, CT, GradientScratch{TT,CT,eltype(plan.rules)})
    try
        sources = PropagationBase._passsources(StorageType(cache), cache, workspace, thread)
        n_tasks = length(sources)
        tasks = PropagationBase._ensurecount!(workspace.tasks, n_tasks)
        for task_id in 1:n_tasks
            scratch = tasks[task_id].scratch
            fill!(resize!(scratch.gradient, n_rotations), 0.0)
            resize!(scratch.rules, n_rotations)
            for k in 1:n_rotations
                scratch.rules[k] = eltype(scratch.rules)(plan.rules[k], scratch.gradient, k)
            end
        end

        _applylayerpass!(cache, plan, gate_of, truncates_each, truncfunc, workspace, sources, thread)

        # the angle of a gate enters each of its rotations, with the factor that the gate gives it
        factors = _anglefactors(layer, thetas, nsites(cache))
        for task_id in 1:n_tasks, k in 1:n_rotations
            gate_grads[gate_of[k]] += factors[k] * tasks[task_id].scratch.gradient[k]
        end
    finally
        PropagationBase._putbackworkspace!(workspace)
    end
    return gate_grads
end

# the sign with which a rotation maps a string to its partner, which a `_GradientRule` pairs the two with
PauliPropagation._partnersign(rule::_RotationRule, term) = _productsign(rule, term)

# The derivative of the angle of every Majorana rotation of the layer by the angle of its gate, in the order in which
# `_classplan` lists the rotations: a fermionic rotation gives each of its rotations twice its coefficient.
function _anglefactors(layer::GateLayer, thetas, n_sites)
    factors = Float64[]
    pushfactors!(gate::MajoranaRotation, theta) = push!(factors, 1.0)
    pushfactors!(gate::FermionicRotation, theta) = append!(factors, 2.0 .* getmajoranarotations(gate, n_sites)[2])
    PropagationBase._foreachgate(pushfactors!, layer, thetas)
    return factors
end


# Gradient contribution of one Majorana rotation: real((i/2) * dual_sum(commutator(generator, op_sum))).
# For a string that anticommutes with the generator, (i/2) times the commutator is its product with the generator,
# with the sign that the rotation gives that product (`_rotationproduct_evengate`). So this is a single pass over
# the strings of the operator, each adding its coefficient times the dual coefficient of its parent string, with
# that sign. The parent is found by looking it up in the merged sum.
# ponytail: one lookup per anticommuting string; port the galloping walk of PauliPropagation's `_commutatordotrange`
# for sorted arrays if this dominates a profile
function _generatorcommutatordot(gate::MajoranaRotation, cache; thread::Bool=true)
    # a lookup needs the sum merged
    merge!(cache; thread)

    gate_ms = _branchmask(gate, cache)
    gate_ms_ps = compute_parity_bits_and_shift(gate_ms, 2 * nfermions(cache))
    return _commutatordot(activesum(cache), cache, gate_ms, gate_ms_ps, omega_L_mult(gate_ms); thread)
end

# in a function of its own so that the closure only captures arguments
function _commutatordot(op_sum, cache, gate_ms, gate_ms_ps, omega_l_gate; thread::Bool=true)
    function commutatoroverlap(term, coeff)
        if _commutes_evengate(term, gate_ms)
            return 0.0
        end
        parent_term, sign = _rotationproduct_evengate(term, gate_ms, gate_ms_ps, omega_l_gate)
        return real(sign * coeff.coeff * getmergedcoeff(op_sum, parent_term).dual)
    end
    return mapreduce(commutatoroverlap, +, cache; init=0.0, thread)
end
