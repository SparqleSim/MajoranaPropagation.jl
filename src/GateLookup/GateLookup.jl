module GateLookup

using MajoranaPropagation
using PauliPropagation.PropagationBase

include("./DataTypes.jl")
include("./create_table.jl")
include("./gates.jl")

export
    MajoranaTransferMap,
    FermionicRotationLookup
end
