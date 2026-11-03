# Monte Carlo

Monte Carlo propagation trades a bounded number of Majorana strings for a statistical error.
`mcpropagate` propagates exactly and resamples the sum down whenever it grows past `max_size`; `mcsample` keeps one random branch per Majorana string and rotation, so the number of strings never changes.
The loops and the resampling strategies are those of `PauliPropagation.PropagationBase`; the Majorana methods add the Heisenberg picture and the sampling rule of a Majorana rotation.

## Majorana methods

```@autodocs
Modules = [MajoranaPropagation]
Pages = ["src/montecarlo.jl"]
```

## Generic interface

```@autodocs
Modules = [MajoranaPropagation.PauliPropagation.PropagationBase]
Pages = ["src/Base/MonteCarlo/mcpropagate.jl", "src/Base/MonteCarlo/mcsample.jl", "src/Base/MonteCarlo/resample.jl"]
```
