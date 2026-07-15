using MajoranaPropagation
using Test
using Random

@testset "MajoranaPropagation.jl" begin
    include("test_commutation_relations.jl")
    include("test_algebra.jl")
    include("check_evengates.jl")
    include("compare_jw.jl")
    include("test_vector.jl")
    include("test_xormerge.jl")
end