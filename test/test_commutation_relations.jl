using Test
using Random
Random.seed!(42)

@testset "Commutation relations" begin
    @testset "single site" begin
        nf = 1
        gamma = MajoranaString(nf, [1])
        gamma_prime = MajoranaString(nf, [2])
        density = MajoranaString(nf, [1, 2])

        @test omega_mult(gamma, gamma_prime) == 1
        @test omega_mult(gamma, density) == 1
    end

    @testset "multi-sites" begin
        nf = 20
        gamma1 = MajoranaString(nf, [1])
        gamma_prime1 = MajoranaString(nf, [2])
        density1 = MajoranaString(nf, [1, 2])
        density2 = MajoranaString(nf, [3, 4])
        gamma2 = MajoranaString(nf, [3])
        gamma_prime2 = MajoranaString(nf, [4])
        gamma3 = MajoranaString(nf, [5])
        gamma_prime3 = MajoranaString(nf, [6])

        #check products
        pref, g1p3 = ms_mult(gamma1, gamma_prime3)
        @test pref == -1im
        pref, p1g2 = ms_mult(gamma_prime1, gamma2)
        @test pref == -1im
        pref, g1p2 = ms_mult(gamma1, gamma_prime2)

        #check single site test in multi sites
        @test omega_mult(gamma1, gamma_prime1) == 1
        @test omega_mult(density1, gamma_prime1) == 1

        @test omega_mult(density1, density2) == 0
        @test omega_mult(gamma1, density2) == 0
        @test omega_mult(gamma_prime1, density2) == 0

        @test omega_mult(g1p3, p1g2) == 0
        @test omega_mult(g1p3, gamma2) == 0
        @test omega_mult(g1p3, gamma1) == 1
        @test omega_mult(g1p3, g1p2) == 1


        ms1 = MajoranaString(nf, [1, 4, 6, 10, 11])
        ms2 = MajoranaString(nf, [2, 15])
        ms3 = MajoranaString(nf, [3, 5, 7, 8, 9])

        @test omega_mult(ms1, ms2) == 0
        @test omega_mult(ms1, ms3) == 1

    end
end