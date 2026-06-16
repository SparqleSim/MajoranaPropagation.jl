module GateLookup

using MajoranaPropagation
using PauliPropagation.PropagationBase
using Combinatorics

include("./DataTypes.jl")
export
    LookupCT
include("./gates_lookup.jl")
include("./create_table.jl")
export
    FermionicRotationLookup
include("./gates.jl")
end