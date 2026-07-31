module FermionToQubitMappings

using MajoranaPropagation
using PauliPropagation

include("./JordanWigner.jl")
export 
    JordanWigner,
    majoranapropagation2yao

end
