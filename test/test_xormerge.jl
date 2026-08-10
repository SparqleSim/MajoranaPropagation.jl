using MajoranaPropagation
using PauliPropagation
using Test
using Random

# Tests for the parallel XOR block-swap tail sort and the gate-aware tail merge
# (src/propagationcache.jl). The parallel chunking only spawns tasks when Julia runs
# with more than one thread; run the suite with JULIA_NUM_THREADS>1 to exercise it.

const xorm_rng = MersenneTwister(42)

# random value of any unsigned width, built from UInt64 chunks (rand(::Type{TT}) is not
# defined for all BitIntegers types)
function xorm_randwide(::Type{TT}) where {TT}
    v = zero(TT)
    for _ in 1:cld(8 * sizeof(TT), 64)
        v = (v << 64) | (rand(xorm_rng, UInt64) % TT)
    end
    return v
end

function xorm_uniqsorted(::Type{TT}, m; mask::TT=~zero(TT)) where {TT}
    s = Set{TT}()
    while length(s) < m
        push!(s, xorm_randwide(TT) & mask)
    end
    return sort!(collect(s))
end

# sort tail = srcs .⊻ g (srcs ascending and unique) with _xorsorttail! and compare
# against Base.sort; also checks coefficients travel with their terms and that the
# result lands in pair A iff popcount(g) is even
function xorm_checksort(srcs::Vector{TT}, g::TT; thread::Bool) where {TT}
    m = length(srcs)
    tail = srcs .⊻ g
    coeffs = collect(Float64, 1:m)
    at = copy(tail)
    ac = copy(coeffs)
    a_terms = view(at, 1:m)
    a_coeffs = view(ac, 1:m)
    b_terms = view(Base.similar(tail), 1:m)
    b_coeffs = view(Base.similar(coeffs), 1:m)
    res_terms, res_coeffs = MajoranaPropagation._xorsorttail!(a_terms, a_coeffs, b_terms, b_coeffs, g; thread)
    p = sortperm(tail)
    return collect(res_terms) == tail[p] &&
           collect(res_coeffs) == coeffs[p] &&
           (iseven(count_ones(g)) ? res_terms === a_terms : res_terms === b_terms)
end

@testset "XOR tail sort vs Base.sort" begin
    for TT in (UInt8, UInt16, UInt64, getinttype(80), getinttype(600))
        nb = 8 * sizeof(TT)
        gs = TT[
            (one(TT) << (nb ÷ 3)) | (one(TT) << (2 * nb ÷ 3)),                        # popcount 2, mid bits
            (one(TT) << (nb - 1)) | one(TT),                                          # popcount 2, top bit (hishift == nbits)
            (one(TT) << (nb - 1)) | (one(TT) << (nb ÷ 2)) | (one(TT) << 3) | one(TT), # popcount 4 incl. top bit
            one(TT),                                                                  # popcount 1 (odd, result in pair B)
            (one(TT) << 5) | (one(TT) << 2) | one(TT),                                # popcount 3 (odd)
        ]
        for g in gs, m in (1, 2, 5, 100, 200)
            srcs = xorm_uniqsorted(TT, m)
            @test xorm_checksort(srcs, g; thread=true)
            @test xorm_checksort(srcs, g; thread=false)
        end
    end

    # sizes crossing _MIN_ELEMS_PER_TASK (= 16384): with >1 thread these run multi-task,
    # so groups straddle chunk boundaries and the binary-search recovery path is exercised
    let TT = UInt64
        g2 = (one(TT) << 40) | (one(TT) << 21)
        g4 = g2 | (one(TT) << 63) | (one(TT) << 3)
        for m in (16_383, 16_385, 40_000, 150_000), g in (g2, g4)
            @test xorm_checksort(xorm_uniqsorted(TT, m), g; thread=true)
        end
        @test xorm_checksort(xorm_uniqsorted(TT, 40_000), g2; thread=false)

        # adversarial structures at multi-task sizes:
        # (i) all sources below 2^18 -> one giant group on every high-bit pass
        for m in (40_000, 150_000), g in (g2, g4)
            srcs = xorm_uniqsorted(TT, m; mask=(one(TT) << 18) - one(TT))
            @test xorm_checksort(srcs, g; thread=true)
        end
        # (ii) dense consecutive sources -> per-element-scale groups on low-bit passes
        @test xorm_checksort(collect(TT, 1:150_000), g2; thread=true)
        @test xorm_checksort(collect(TT, 1:150_000), TT(0b11); thread=true)
    end
end

# build a merge-ready cache: sorted+deduped head of length n_head, then the tail appended
# behind it, with the sorted prefix set to n_head (exactly the state _mergeafterapply! sees)
function xorm_buildcache(head_terms::Vector{TT}, head_coeffs, tail_terms, tail_coeffs) where {TT}
    terms = vcat(head_terms, tail_terms)
    coeffs = vcat(head_coeffs, tail_coeffs)
    vms = VectorMajoranaSum(4 * sizeof(TT), false, terms, coeffs, length(head_terms))
    return MajoranaPropagation.VectorMajoranaPropagationCache(vms)
end

# reference: dict accumulation of head + tail, sorted by term, optionally dropping
# collided keys per truncfunc (matching the merge kernel's collision-only truncation)
function xorm_refmerge(head_terms::Vector{TT}, head_coeffs, tail_terms, tail_coeffs; truncfunc=nothing) where {TT}
    d = Dict{TT,Float64}()
    collided = Set{TT}()
    for (t, c) in zip(head_terms, head_coeffs)
        d[t] = get(d, t, 0.0) + c
    end
    for (t, c) in zip(tail_terms, tail_coeffs)
        haskey(d, t) && push!(collided, t)
        d[t] = get(d, t, 0.0) + c
    end
    if !isnothing(truncfunc)
        for t in collided
            truncfunc(t, d[t]) && delete!(d, t)
        end
    end
    ks = sort!(collect(keys(d)))
    return ks, [d[k] for k in ks]
end

function xorm_checkmerge(cache, ref_terms, ref_coeffs)
    n = MajoranaPropagation.activesize(cache)
    got_terms = collect(MajoranaPropagation.activeterms(cache))
    got_coeffs = collect(MajoranaPropagation.activecoeffs(cache))
    return n == length(ref_terms) &&
           got_terms == ref_terms &&
           got_coeffs == ref_coeffs &&
           MajoranaPropagation.sortedprefix(MajoranaPropagation.mainsum(cache)) == n
end

@testset "xorsortedtailmerge! unit" begin
    TT = UInt16
    g = (one(TT) << 9) | (one(TT) << 2)   # popcount 2 (even, physical case)
    g_odd = one(TT) << 4                   # popcount 1 (odd: sorted tail lands in pair B)

    # head containing pairs (a, a ⊻ g) so that rotating a collides with an existing term
    base = xorm_uniqsorted(TT, 40)
    head_terms = sort!(unique!(vcat(base, base[1:10] .⊻ g)))
    head_coeffs = randn(xorm_rng, length(head_terms))
    n_head = length(head_terms)

    # ascending subset of the head as rotation sources (mirrors _applymajoranarotation!)
    sources = head_terms[1:2:end]
    tail_coeffs = randn(xorm_rng, length(sources))

    for (glabel, gg) in (("even", g), ("odd", g_odd))
        tail_terms = sources .⊻ gg

        # fresh-alloc branch: the constructor leaves aux with zero headroom
        cache = xorm_buildcache(head_terms, head_coeffs, tail_terms, tail_coeffs)
        MajoranaPropagation.xorsortedtailmerge!(cache, gg)
        ref_t, ref_c = xorm_refmerge(head_terms, head_coeffs, tail_terms, tail_coeffs)
        @test xorm_checkmerge(cache, ref_t, ref_c)

        # aux-headroom branch: grow the cache so aux has spare capacity beyond n_new
        cache = xorm_buildcache(head_terms, head_coeffs, tail_terms, tail_coeffs)
        resize!(cache, 4 * (n_head + length(tail_terms)))
        MajoranaPropagation.xorsortedtailmerge!(cache, gg)
        @test xorm_checkmerge(cache, ref_t, ref_c)
    end

    # truncfunc fires only on actual collisions: engineer an exact cancellation
    let tail_terms = sources .⊻ g
        cancel_term = tail_terms[3]                      # == sources[3] ⊻ g
        ht = sort!(unique!(vcat(head_terms, cancel_term)))
        hc = randn(xorm_rng, length(ht))
        hc[searchsortedfirst(ht, cancel_term)] = 0.5
        tc = copy(tail_coeffs)
        tc[3] = -0.5
        truncfunc = (trm, coeff) -> abs(coeff) < 1e-12

        cache = xorm_buildcache(ht, hc, tail_terms, tc)
        MajoranaPropagation.xorsortedtailmerge!(cache, g; truncfunc)
        ref_t, ref_c = xorm_refmerge(ht, hc, tail_terms, tc; truncfunc)
        @test cancel_term ∉ ref_t
        @test xorm_checkmerge(cache, ref_t, ref_c)
    end

    # fallback routes must produce the same result via the generic merge!
    let tail_terms = sources .⊻ g
        ref_t, ref_c = xorm_refmerge(head_terms, head_coeffs, tail_terms, tail_coeffs)
        cache = xorm_buildcache(head_terms, head_coeffs, tail_terms, tail_coeffs)
        MajoranaPropagation.xorsortedtailmerge!(cache, nothing)          # no gate string
        @test xorm_checkmerge(cache, ref_t, ref_c)
        cache = xorm_buildcache(head_terms, head_coeffs, tail_terms, tail_coeffs)
        MajoranaPropagation.xorsortedtailmerge!(cache, UInt64(g))        # eltype mismatch
        @test xorm_checkmerge(cache, ref_t, ref_c)

        # gate weight > 4: XOR preconditions hold but the weight cutoff forces the fallback
        g6 = (one(TT) << 11) | (one(TT) << 9) | (one(TT) << 7) | (one(TT) << 5) | (one(TT) << 2) | one(TT)
        @test get_weight(g6) == 6
        tail6 = sources .⊻ g6
        ref_t6, ref_c6 = xorm_refmerge(head_terms, head_coeffs, tail6, tail_coeffs)
        cache = xorm_buildcache(head_terms, head_coeffs, tail6, tail_coeffs)
        MajoranaPropagation.xorsortedtailmerge!(cache, g6)
        @test xorm_checkmerge(cache, ref_t6, ref_c6)
    end
end

# O(N) agreement check (the AbstractTermSum `==` re-merges copies and looks terms up
# one by one, which is far too slow at the term counts these blocks reach)
function xorm_sums_agree(msum, vmsum; atol=1e-10)
    d_dict = Dict(k => v for (k, v) in msum)
    d_vec = Dict(k => v for (k, v) in vmsum)
    return length(d_dict) == length(d_vec) &&
           all(isapprox(v, get(d_vec, k, Inf); atol) for (k, v) in d_dict)
end

@testset "XOR merge end-to-end (dict vs vector)" begin
    # interacting 2D problem large enough that per-rotation tails exceed
    # 2 * _MIN_ELEMS_PER_TASK, so threaded runs merge with multiple tasks
    nx, ny = 5, 4
    topo = rectangletopology(nx, ny)
    nf = nx * ny
    circ = FermionicRotation[]
    thetas = Float64[]
    for pair in topo
        push!(circ, FermionicRotation(:hop, pair))
        push!(thetas, 0.1)
        push!(circ, FermionicRotation(:nn, pair))
        push!(thetas, 0.1)
    end

    obs = MajoranaSum(nf, :n, (nf + 1) ÷ 2)
    obs_vec = VectorMajoranaSum(deepcopy(obs))
    for _ in 1:4
        propagate!(circ, obs, thetas; min_abs_coeff=1e-6)
        propagate!(circ, obs_vec, thetas; min_abs_coeff=1e-6)
        @test length(obs) == length(obs_vec)
        @test xorm_sums_agree(obs, obs_vec)
    end
    @test length(obs_vec) > 100_000

    # bare MajoranaRotation circuit: exercises the applymergetruncate! overload that
    # routes bare rotations through the gate-aware merge
    nf = 10
    TT = getinttype(nf)
    mr_circ = MajoranaRotation{TT}[]
    mr_thetas = Float64[]
    for _ in 1:40
        w = rand(xorm_rng, (2, 4))
        gammas = Int[]
        while length(gammas) < w
            gg = rand(xorm_rng, 1:2*nf)
            gg in gammas || push!(gammas, gg)
        end
        push!(mr_circ, MajoranaRotation(MajoranaString(nf, gammas)))
        push!(mr_thetas, randn(xorm_rng))
    end
    obs = MajoranaSum(nf, :n, 5)
    obs_vec = VectorMajoranaSum(deepcopy(obs))
    for _ in 1:2
        propagate!(mr_circ, obs, mr_thetas; min_abs_coeff=1e-6)
        propagate!(mr_circ, obs_vec, mr_thetas; min_abs_coeff=1e-6)
        @test length(obs) == length(obs_vec)
        @test xorm_sums_agree(obs, obs_vec)
    end
end
