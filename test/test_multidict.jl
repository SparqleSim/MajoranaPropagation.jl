using Test
using MajoranaPropagation
using PauliPropagation

@testset "MajoranaMultiPropagationCache" begin
    seed = MajoranaSum(Float64, 2, [1], false)
    key = only(keys(seed.Majoranas))

    multi = MajoranaSumMulti(seed, _ -> 1, 1)
    cache = MajoranaMultiPropagationCache(multi)

    main = mainsum(cache)
    aux = auxsum(cache)
    empty!(main.MultiMajoranas[1])
    empty!(aux.MultiMajoranas[1])

    main.MultiMajoranas[1][key] = 2.0
    aux.MultiMajoranas[1][key] = 3.0

    merge!(cache; level_mapper = _ -> 1)

    @test length(cache) == 1
    @test main.MultiMajoranas[1][key] == 5.0
    @test isempty(aux.MultiMajoranas[1])
end