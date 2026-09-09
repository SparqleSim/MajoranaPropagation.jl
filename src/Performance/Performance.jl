module Performance
# Opt-in module with fused versions of `applymergetruncate!` that truncate while a gate is applied
# instead of in a separate `truncate!` pass, mirroring `PauliPropagation.Performance`.
# Opt in via `Performance.propagate` / `Performance.propagate!`, or by passing `fused=true` to the
# regular propagation functions.

using MajoranaPropagation
using PauliPropagation
using PauliPropagation.PropagationBase

using AcceleratedKernels
const AK = AcceleratedKernels

const MP = MajoranaPropagation

# fused overloads for `VectorMajoranaSum` caches and the fermionic gate wrappers
include("./fused_vector.jl")

# fused overload for the dictionary-based `MajoranaSum` cache
include("./fused_dict.jl")

"""
    propagate(circuit, thing, thetas=nothing; fused::Bool=true, kwargs...)

Like `MajoranaPropagation.propagate`, but defaults `fused=true` to use this module's fused `applymergetruncate!` overloads.
Pass `fused=false` for the default behavior.
The fused kernels assume that `thing` already satisfies the truncations, as it does after any `truncate!` or `propagate!` call; call `truncate!` first otherwise.
"""
function propagate(circuit, thing, thetas=nothing; fused::Bool=true, kwargs...)
    return MP.propagate(circuit, thing, thetas; fused, kwargs...)
end

"""
    propagate!(circuit, thing, thetas=nothing; fused::Bool=true, kwargs...)

In-place counterpart of `propagate`. See `propagate` for details.
"""
function propagate!(circuit, thing, thetas=nothing; fused::Bool=true, kwargs...)
    return MP.propagate!(circuit, thing, thetas; fused, kwargs...)
end

end
