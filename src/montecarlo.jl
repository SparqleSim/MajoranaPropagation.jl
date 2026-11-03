###
##
# Monte Carlo propagation of Majorana sums. `mcpropagate!` runs the deterministic gates and resamples
# the sum down whenever it grows past `max_size`; `mcsample!` keeps one random branch per Majorana
# string and rotation. The loops and the resampling strategies are PropagationBase's; this file
# converts the circuit to the Heisenberg picture and writes the sampling rule of a Majorana rotation.
##
###

"""
    mcpropagate!(circuit, prop_cache::AbstractMajoranaPropagationCache, thetas=nothing; max_size, heisenberg=true, thread=true, kwargs...)

In-place Monte Carlo propagation of a Majorana propagation cache through `circuit` (see `mcpropagate`), by default in the Heisenberg picture.
`kwargs` are passed to `applymergetruncate!` (e.g. `min_abs_coeff`, `max_weight`) and to the resampling strategy (e.g. `squared`).
Coefficients wrapped in a `MajoranaFrequencyTracker` cannot be resampled, and resampling may drop the identity that `normalize_coeffs=true` divides by.
"""
function PropagationBase.mcpropagate!(circuit, prop_cache::AbstractMajoranaPropagationCache, thetas=nothing; heisenberg=true, kwargs...)
    circuit, thetas = _preparecircuit(circuit, thetas, heisenberg)
    return PropagationBase._propagate!(PropagationBase.applymergetruncateresample!, circuit, prop_cache, thetas; kwargs...)
end

"""
    mcsample!(circuit, msum::AbstractMajoranaSum, thetas=nothing; squared=false, heisenberg=true, thread=true, kwargs...)

In-place Monte Carlo path sampling of a Majorana sum through `circuit` (see `mcsample`), by default in the Heisenberg picture.
Every Majorana string keeps one random branch of each rotation, so the number of strings never changes.
A `VectorMajoranaSum` keeps walkers that land on the same string apart, a `MajoranaSum` merges them.
Imaginary-time rotations are not supported.
Sample a `MultiMajoranaSum` itself, not its propagation cache.
"""
function PropagationBase.mcsample!(circuit, msum::AbstractMajoranaSum, thetas=nothing; heisenberg=true, kwargs...)
    circuit, thetas = _preparecircuit(circuit, thetas, heisenberg)
    return PropagationBase._propagate!(PropagationBase.mcapplytoall!, circuit, msum, thetas; kwargs...)
end

# PropagationBase's cache method rebuilds the zones as views, which cannot be resized
function PropagationBase.mcsample!(circuit, prop_cache::MultiMajoranaPropagationCache, thetas=nothing; kwargs...)
    throw(ArgumentError("mcsample! is not supported on a MultiMajoranaPropagationCache; sample the MultiMajoranaSum itself"))
end

"""
    mcapplytoall!(gate::MajoranaRotation, msum::AbstractMajoranaSum, theta; squared=false, thread=true, kwargs...)
    mcapplytoall!(gate::FermionicRotation, msum::AbstractMajoranaSum, theta; squared=false, thread=true, kwargs...)

Stochastically apply one rotation to every Majorana string of `msum` in place, the Monte Carlo analogue of `applytoall!`.
A string that anticommutes with the rotation keeps its sine branch with probability `|sin|/(|sin|+|cos|)` and its cosine branch otherwise, reweighted by `|sin|+|cos|` so that the sum stays unbiased in expectation.
With `squared=true` the probabilities are `sin^2` and `cos^2`, which preserves `|coeff|`.
A `MajoranaFrequencyTracker` counts the branch a string keeps.
A `FermionicRotation` samples each of its Majorana rotations in turn.
"""
function PropagationBase.mcapplytoall!(gate::MajoranaRotation, msum::AbstractMajoranaSum, theta; squared::Bool=false, thread::Bool=true, kwargs...)
    gate_ms = _branchmask(gate, msum)
    gate_ms_ps = compute_parity_bits_and_shift(gate_ms, 2 * nfermions(msum))
    omega_l_gate = omega_L_mult(gate_ms)

    power = squared ? 2 : 1
    sin_val, cos_val = sin(theta), cos(theta)
    # a branch is kept with probability proportional to |sin|^power or |cos|^power, and carries the
    # normalization so that the sum stays unbiased
    normalization = abs(sin_val)^power + abs(cos_val)^power
    p_sin = abs(sin_val)^power / normalization
    sin_weight = normalization * sign(sin_val)^power
    cos_weight = normalization * sign(cos_val)^power

    function sample(term, coeff)
        _commutes_evengate(term, gate_ms) && return (term, coeff)
        rand() < p_sin || return (term, _applycos(coeff, cos_weight))
        new_term, prod_sign = _rotationproduct_evengate(term, gate_ms, gate_ms_ps, omega_l_gate)
        return (new_term, _applysin(coeff, sin_weight * prod_sign^power))
    end

    return map!(sample, msum; thread)
end

function PropagationBase.mcapplytoall!(gate::FermionicRotation, msum::AbstractMajoranaSum, theta; kwargs...)
    ms_rotations, coeffs, _ = getmajoranarotations(gate, nsites(msum))
    for (gate_ms, coeff) in zip(ms_rotations, coeffs)
        # `MajoranaRotation` implements exp(-i * theta/2 * mstring), hence the factor 2
        mcapplytoall!(gate_ms, msum, theta * coeff * 2.0; kwargs...)
    end
    return msum
end
