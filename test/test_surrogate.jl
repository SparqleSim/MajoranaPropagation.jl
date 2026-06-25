# The surrogate builds the symbolic graph once and is evaluated for many parameter vectors.
# Ground truth is the ordinary numeric Majorana propagation at the same parameters.
@testset "Surrogate" begin

    @testset "MajoranaRotation circuit" begin
        n = 4                      # spinless: n sites = n fermions, 2n Majorana modes
        obs_site = 2
        fock = FockState(n, [1, 3])

        # circuit of plain MajoranaRotations (even-weight Majorana strings)
        circ = MajoranaRotation[
            MajoranaRotation(MajoranaString(n, [1, 4])),
            MajoranaRotation(MajoranaString(n, [3, 6])),
            MajoranaRotation(MajoranaString(n, [2, 5])),
            MajoranaRotation(MajoranaString(n, [4, 7])),
            MajoranaRotation(MajoranaString(n, [1, 2])),
        ]
        nparams = countparameters(circ)

        # build the surrogate once
        surr = propagate(circ, MajoranaPropagation.wrapcoefficients(MajoranaSum(n, :n, obs_site), MajoranaNodePathProperties))

        Random.seed!(1)
        for _ in 1:6
            thetas = randn(nparams)
            ref = overlapwithfock(propagate(circ, MajoranaSum(n, :n, obs_site), thetas; min_abs_coeff=-1.0), fock)
            evaluate!(surr, thetas)               # reuses the same surrogate across parameter sets
            @test isapprox(overlapwithfock(surr, fock), ref; atol=1e-10)
        end
    end

    @testset "FermionicRotation circuit" begin
        n = 6
        obs_site = 3
        topo = bricklayertopology(n)
        fock = FockState(n, collect(1:2:n))

        circ = FermionicRotation[]
        for (i, j) in topo
            push!(circ, FermionicRotation(:hop, [i, j]))
        end
        for (i, j) in topo
            push!(circ, FermionicRotation(:nn, [i, j]))
        end
        nparams = countparameters(circ)

        surr = propagate(circ, MajoranaPropagation.wrapcoefficients(MajoranaSum(n, :n, obs_site), MajoranaNodePathProperties))

        Random.seed!(2)
        for _ in 1:4
            thetas = randn(nparams)
            ref = overlapwithfock(propagate(circ, MajoranaSum(n, :n, obs_site), thetas; min_abs_coeff=-1.0), fock)
            evaluate!(surr, thetas)
            @test isapprox(overlapwithfock(surr, fock), ref; atol=1e-10)
        end
    end

    @testset "max_weight truncation matches numeric" begin
        n = 6
        obs_site = 3
        max_weight = 4
        topo = bricklayertopology(n)
        fock = FockState(n, collect(1:2:n))

        circ = FermionicRotation[]
        for (i, j) in topo
            push!(circ, FermionicRotation(:hop, [i, j]))
        end
        nparams = countparameters(circ)

        surr = propagate(circ, MajoranaPropagation.wrapcoefficients(MajoranaSum(n, :n, obs_site), MajoranaNodePathProperties); max_weight=max_weight)

        Random.seed!(3)
        for _ in 1:4
            thetas = randn(nparams)
            ref = overlapwithfock(propagate(circ, MajoranaSum(n, :n, obs_site), thetas; min_abs_coeff=-1.0, max_weight=max_weight), fock)
            evaluate!(surr, thetas)
            @test isapprox(overlapwithfock(surr, fock), ref; atol=1e-10)
        end
    end
end
