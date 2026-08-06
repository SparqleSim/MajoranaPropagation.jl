using MajoranaPropagation
using PauliPropagation
using Test
using Random

# MajoranaFrequencyTracker on both backends. The vector sums are constructed with a full
# sorted prefix so every merge takes the counter-exact sorted-tail path (the generic
# fallback merge resets counters, see the wrapcoefficients docstring).

const ft_rng = MersenneTwister(7)

ft_agree(a, b; atol=1e-10) =
    isapprox(a.coeff, b.coeff; atol) && a.freq == b.freq && a.nsins == b.nsins && a.ncos == b.ncos

function ft_sums_agree(msum, vmsum; atol=1e-10)
    d1 = Dict(k => v for (k, v) in msum)
    d2 = Dict(k => v for (k, v) in vmsum)
    return length(d1) == length(d2) &&
           all(haskey(d2, k) && ft_agree(v, d2[k]; atol) for (k, v) in d1)
end

# sorted vector copy of a dict sum with the sorted prefix covering all terms
function ft_sorted_vector(msum)
    ks = sort!(collect(keys(msum.Majoranas)))
    cs = [msum.Majoranas[k] for k in ks]
    return VectorMajoranaSum(msum.nsites, msum.is_spinful, ks, cs, length(ks))
end

@testset "MajoranaFrequencyTracker" begin
    # constructors, zero, numcoefftype
    @test zero(MajoranaFrequencyTracker{Float64}) === MajoranaFrequencyTracker(0.0, 0, 0, 0)
    @test zero(MajoranaFrequencyTracker{ComplexF64}).coeff === 0.0 + 0.0im
    @test MajoranaFrequencyTracker(1) === MajoranaFrequencyTracker(1.0, 0, 0, 0)
    @test numcoefftype(MajoranaFrequencyTracker{Float64}) == Float64

    nf = 4
    obs = MajoranaSum(nf, :n, 2)

    # wrap/unwrap round trip, dict and vector
    wobs = wrapcoefficients(obs, MajoranaFrequencyTracker)
    @test all(v === MajoranaFrequencyTracker(Float64(v.coeff)) for (k, v) in wobs)
    @test Dict(k => v for (k, v) in unwrapcoefficients(wobs)) == Dict(k => v for (k, v) in obs)

    wvobs = wrapcoefficients(ft_sorted_vector(obs), MajoranaFrequencyTracker)
    @test wvobs._terms_sorted == length(obs)
    @test unwrapcoefficients(wvobs)._terms_sorted == length(obs)
    @test Dict(k => v for (k, v) in unwrapcoefficients(wvobs)) == Dict(k => v for (k, v) in obs)

    # hand-computed single rotation: gate g1g2 anticommutes with the term g2g3
    theta = 0.3
    term_old = MajoranaString(nf, [2, 3]).gammas
    term_new = term_old ⊻ MajoranaString(nf, [1, 2]).gammas
    gate = MajoranaRotation(MajoranaString(nf, [1, 2]))

    wsum = wrapcoefficients(MajoranaSum(nf, false, Dict(term_old => 1.0)), MajoranaFrequencyTracker)
    propagate!([gate], wsum, [theta])
    d = Dict(k => v for (k, v) in wsum)
    @test length(d) == 2
    @test ft_agree(d[term_old], MajoranaFrequencyTracker(cos(theta), 1, 0, 1))
    @test abs(d[term_new].coeff) ≈ abs(sin(theta))
    @test (d[term_new].freq, d[term_new].nsins, d[term_new].ncos) == (1, 1, 0)

    wvsum = wrapcoefficients(ft_sorted_vector(MajoranaSum(nf, false, Dict(term_old => 1.0))), MajoranaFrequencyTracker)
    propagate!([gate], wvsum, [theta])
    @test ft_sums_agree(wsum, wvsum)

    # backend agreement on a random even-weight <= 4 rotation circuit
    nf = 6
    TT = getinttype(nf)
    circ = MajoranaRotation{TT}[]
    thetas = Float64[]
    for _ in 1:30
        w = rand(ft_rng, (2, 4))
        gammas = Int[]
        while length(gammas) < w
            gg = rand(ft_rng, 1:2*nf)
            gg in gammas || push!(gammas, gg)
        end
        push!(circ, MajoranaRotation(MajoranaString(nf, gammas)))
        push!(thetas, randn(ft_rng))
    end
    obs = MajoranaSum(nf, :n, 3)

    wobs = wrapcoefficients(obs, MajoranaFrequencyTracker)
    wvobs = wrapcoefficients(ft_sorted_vector(obs), MajoranaFrequencyTracker)
    propagate!(circ, wobs, thetas; min_abs_coeff=1e-8)
    propagate!(circ, wvobs, thetas; min_abs_coeff=1e-8)
    @test length(wobs) > 2
    @test ft_sums_agree(wobs, wvobs)
    @test any(v.freq > 0 for (k, v) in wvobs)

    # max_freq/max_sins truncation fires on both backends and they agree
    tobs = wrapcoefficients(obs, MajoranaFrequencyTracker)
    tvobs = wrapcoefficients(ft_sorted_vector(obs), MajoranaFrequencyTracker)
    propagate!(circ, tobs, thetas; min_abs_coeff=1e-8, max_freq=2, max_sins=1)
    propagate!(circ, tvobs, thetas; min_abs_coeff=1e-8, max_freq=2, max_sins=1)
    @test ft_sums_agree(tobs, tvobs)
    @test length(tvobs) < length(wvobs)
    @test all(v.freq <= 2 && v.nsins <= 1 for (k, v) in tvobs)

    # reset_tracker! on both backends
    reset_tracker!(wobs)
    reset_tracker!(wvobs)
    @test all(v.freq == 0 && v.nsins == 0 && v.ncos == 0 for (k, v) in wobs)
    @test all(v.freq == 0 && v.nsins == 0 && v.ncos == 0 for (k, v) in wvobs)
    @test ft_sums_agree(wobs, wvobs)
end
