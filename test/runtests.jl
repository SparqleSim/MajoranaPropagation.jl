using MajoranaPropagation
using PauliPropagation
using Test
using Random

@testset "MajoranaPropagation.jl" begin
    include("test_commutation_relations.jl")
    include("test_algebra.jl")
    include("compare_jw.jl")
end