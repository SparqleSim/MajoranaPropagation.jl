using MajoranaPropagation
using Test
using Random

@testset "MajoranaPropagation.jl" begin
    # dense Jordan-Wigner reference implementation
    include("testutils_dense.jl")

    include("test_commutation_relations.jl")
    include("test_algebra.jl")
    include("check_evengates.jl")
    include("compare_jw.jl")
    include("test_vector.jl")
    include("test_fock_dense.jl")
    include("test_imaginary.jl")
    include("test_truncations.jl")
    include("test_sum_algebra.jl")
    include("test_constructors.jl")
    include("test_jordanwigner.jl")
    include("test_xormerge.jl")
    include("test_performance.jl")
end