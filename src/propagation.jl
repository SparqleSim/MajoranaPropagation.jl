"""
    propagate(circuit, msum::AbstractMajoranaSum, thetas=nothing; min_abs_coeff=1e-10, max_weight=Inf, max_unpaired=Inf,  unpaired_mask=nothing, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...)
    propagate!(circuit, msum::AbstractMajoranaSum, thetas=nothing; min_abs_coeff=1e-10, max_weight=Inf, max_unpaired=Inf,  unpaired_mask=nothing, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...)

Propagate a Majorana sum `msum` through the circuit `circ`. 
By default this is done in the Heisenberg picture, as indicated by `heisenberg=true`. 
This means that the circuit is applied to the Majorana sum in reverse order, and the action of each gate is its conjugate action.
In `propagate()` the Majorana sum `msum` is deepcopied and passed into the in-place propagation function `propagate!()`.
Parameters for the parametrized gates in `circ` are given by `thetas`, and need to be passed as if the circuit was applied as written in the Schrödinger picture.
If thetas are not passed, the circuit must contain only non-parametrized `StaticGates`.
Default truncations are `min_abs_coeff`, `max_weight`, `max_freq`, and `max_sins`.
`max_freq`, and `max_sins` will lead to automatic conversion if the coefficients are not already wrapped in suitable `PathProperties` objects.
A custom truncation function can be passed as `customtruncfunc` with the signature customtruncfunc(pstr::PauliStringType, coefficient)::Bool.
Further `kwargs` are passed to the lower-level functions `applymergetruncate!`, `applytoall!`, and `apply`.
"""
function PropagationBase.propagate(circuit, msum::AbstractMajoranaSum, thetas=nothing; max_weight=Inf, min_abs_coeff=1e-10, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...)
    CT = coefftype(msum)

    # if max_freq and max_sins are used, and no PathProperties used, automatically wrap the coefficients in `PauliFreqTracker` 
    msum = _check_wrapping_into_paulifreqtracker(msum, max_freq, max_sins)

    # check that max_freq and max_sins are only used a PathProperties type tracking them
    _checkfreqandsinfields(msum, max_freq, max_sins)

    # run the in-place propagation function on a deepcopy of the input psum
    msum = propagate!(circuit, deepcopy(msum), thetas; max_weight, min_abs_coeff, max_freq, max_sins, customtruncfunc, heisenberg, kwargs...)

    # if the input psum was not a `PauliFreqTracker`, and the corresponding truncations were set,we need to unwrap the coefficients
    msum = _check_unwrap_from_paulifreqtracker(CT, msum)

    return msum
end


"""
    propagate!(circuit, prop_cache::AbstractMajoranaPropagationCache, thetas=nothing; min_abs_coeff=1e-10, max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...)

In-place propagation of an `AbstractMajoranaPropagationCache` through the circuit `circ` in the Heisenberg picture.
"""
function PropagationBase.propagate!(circuit, prop_cache::AbstractMajoranaPropagationCache, thetas=nothing; max_weight=Inf, min_abs_coeff=1e-10, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...)

    # if circuit is actually a single gate, promote it to a list [gate]
    # similarly the thetas if it is a single number
    circuit, thetas = PropagationBase._promotecircandparams(circuit, thetas)

    # if thetas is nothing, the circuit must contain only StaticGates
    # also check if the length of thetas equals the number of parametrized gates
    PropagationBase._checknumberofparams(circuit, thetas)

    if heisenberg
        # this usually just reverses circuit and parameter order
        circuit, thetas = toheisenberg(circuit, thetas)
    else
        # this usually entails a conversion of how gates act
        circuit, thetas = toschrodinger(circuit, thetas)
    end

    return PropagationBase._propagate!(circuit, prop_cache, thetas; max_weight, min_abs_coeff, max_freq, max_sins, customtruncfunc, kwargs...)
end