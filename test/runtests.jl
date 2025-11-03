using Pkg
Pkg.activate("./MajoranaPropagation")  # Activate the project environment
#Pkg.activate("..")  # Activate the current directory as the project environment

using Revise
using MajoranaPropagation
using Test
using Random

@testset "MajoranaPropagation.jl" begin
    include("test_commutation_relations.jl")
end