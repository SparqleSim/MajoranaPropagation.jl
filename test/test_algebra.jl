
using Test
using Random

Random.seed!(42)


function compute_parity_bits_and_shift_legacy_version(u::TT, Nbits::Int) where {TT<:Integer}
    p::TT = 0
    parity::TT = 0
    for k in TT(0):TT(Nbits - 2)
        if (u >> k) & 1 == 1
            parity ⊻= 1
        end
        p |= parity << (k + 1)
    end
    return p
end

@testset "Majorana algebra" begin
    @testset "parites" begin
    n_fermions = 10 
    Nbits = 2 * n_fermions

    weights = [2, 4, 8]
    n_tests = 20
    for weight in weights 
        for _=1:n_tests
            sites::Vector{Int64} = []
            while length(sites) < weight
                rand_int = rand(1:Nbits)
                if !(rand_int in sites)
                    push!(sites, rand_int)
                end
            end
            sites = sort(sites)
            ms = MajoranaString(n_fermions, sites)
            r1 = compute_parity_bits_and_shift_legacy_version(ms.gammas, Nbits)
            r2 = MajoranaPropagation.compute_parity_bits_and_shift(ms.gammas, Nbits)
            @test r1 == r2
        end
    end

    n_fermions = 100 
    Nbits = 2 * n_fermions

    weights = [2, 4, 8, 16, 32]
    n_tests = 20
    for weight in weights 
        for _=1:n_tests
            sites::Vector{Int64} = []
            while length(sites) < weight
                rand_int = rand(1:Nbits)
                if !(rand_int in sites)
                    push!(sites, rand_int)
                end
            end
            sites = sort(sites)
            ms = MajoranaString(n_fermions, sites)
            r1 = compute_parity_bits_and_shift_legacy_version(ms.gammas, Nbits)
            r2 = MajoranaPropagation.compute_parity_bits_and_shift(ms.gammas, Nbits)
            @test r1 == r2
        end
    end
end

end