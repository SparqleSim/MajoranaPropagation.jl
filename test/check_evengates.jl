using MajoranaPropagation
using PauliPropagation
using Test
using Random

Random.seed!(42)

# `applytoall!`, `mainsum`, `auxsum` are exported by PauliPropagation (MajoranaPropagation adds
# its methods to the same PropagationBase generics), so `using PauliPropagation` makes them
# reachable unqualified. `commutes` is exported by BOTH packages and must stay qualified.

# random Majorana string integer with the given weight (0 allowed), any parity
function random_gammas(nfermions::Int, weight::Int)
    weight == 0 && return getinttype(nfermions)(0)
    sites = sort!(shuffle!(collect(1:2*nfermions))[1:weight])
    return MajoranaString(nfermions, sites).gammas
end

# random even-weight, nonzero gate string
function random_even_gate(nfermions::Int)
    return random_gammas(nfermions, 2 * rand(1:min(4, nfermions)))
end

@testset "Even-gate kernels" begin

    @testset "even-weight constructor invariant" begin
        nf = 5
        for R in (MajoranaRotation, ImaginaryMajoranaRotation)
            @test_throws ArgumentError R(MajoranaString(nf, [1]))
            @test_throws ArgumentError R(MajoranaString(nf, [2, 3, 7]).gammas)
            @test R(MajoranaString(nf, [1, 2])) isa R
            @test R(MajoranaString(nf, [1, 4, 5, 9]).gammas) isa R
        end
    end

    # the specialized kernels assume an even-weight GATE but must be exact for terms of ANY
    # parity; checked against the general algebra for native and BitIntegers key widths
    @testset "kernels match general algebra (nfermions=$nf)" for nf in (4, 16, 28)
        nbits = 2 * nf
        for trial in 1:500
            gate = random_even_gate(nf)
            gate_ps = MajoranaPropagation.compute_parity_bits_and_shift(gate, nbits)
            omega_l_gate = omega_L_mult(gate)
            # weight 0 (identity) and odd weights are deliberately included
            term = random_gammas(nf, trial <= 5 ? trial - 1 : rand(0:nbits))

            @test MajoranaPropagation._commutes_evengate(term, gate) == MajoranaPropagation.commutes(gate, term)

            if !MajoranaPropagation.commutes(gate, term)
                nt_ref, sg_ref = majoranarotationproduct(term, gate, gate_ps)
                nt, sg = MajoranaPropagation._rotationproduct_evengate(term, gate, gate_ps, omega_l_gate)
                @test nt == nt_ref && sg == sg_ref
            else
                # commuting branch (imaginary time): gate*term == term*gate, prefactor real
                pref_ref, nt_ref = ms_mult(gate, term, nf)
                nt, sg = MajoranaPropagation._rotationproduct_evengate_commuting(term, gate, gate_ps, omega_l_gate)
                @test pref_ref isa Real
                @test nt == nt_ref && sg == pref_ref
            end
        end
    end

    @testset "real-time applytoall! matches general-kernel reference" begin
        nf = 9
        nbits = 2 * nf
        msum = MajoranaSum(nf)
        for _ in 1:300
            set!(msum, random_gammas(nf, rand(0:nbits)), randn())
        end
        theta = 0.7
        for _ in 1:10
            gate_int = random_even_gate(nf)
            gate_ps = MajoranaPropagation.compute_parity_bits_and_shift(gate_int, nbits)

            cache = MajoranaPropagationCache(deepcopy(msum))
            applytoall!(MajoranaRotation(gate_int), cache, theta)
            new_main = Dict(k => v for (k, v) in mainsum(cache))
            new_aux = Dict(k => v for (k, v) in auxsum(cache))

            ref_main = Dict(k => v for (k, v) in msum)
            ref_aux = Dict{keytype(ref_main),Float64}()
            for (k, c) in msum
                MajoranaPropagation.commutes(gate_int, k) && continue
                nt, sg = majoranarotationproduct(k, gate_int, gate_ps)
                ref_main[k] = c * cos(theta)
                ref_aux[nt] = c * sin(theta) * sg
            end
            @test new_main == ref_main
            @test new_aux == ref_aux
        end
    end

    @testset "imaginary-time applytoall! matches ms_mult reference" begin
        nf = 9
        nbits = 2 * nf
        msum = MajoranaSum(nf)
        for _ in 1:300
            set!(msum, random_gammas(nf, rand(0:nbits)), randn())
        end
        beta = 0.37
        for _ in 1:10
            gate_int = random_even_gate(nf)

            cache = MajoranaPropagationCache(deepcopy(msum))
            applytoall!(ImaginaryMajoranaRotation(gate_int), cache, beta)
            new_main = Dict(k => v for (k, v) in mainsum(cache))
            new_aux = Dict(k => v for (k, v) in auxsum(cache))

            ref_main = Dict(k => v for (k, v) in msum)
            ref_aux = Dict{keytype(ref_main),Float64}()
            for (k, c) in msum
                MajoranaPropagation.commutes(gate_int, k) || continue
                pref, nt = ms_mult(gate_int, k, nf)
                ref_main[k] = c * cosh(beta)
                ref_aux[nt] = c * (-sinh(beta)) * real(pref)
            end
            @test new_main == ref_main
            @test new_aux == ref_aux
        end
    end

    @testset "dict and vector backends agree" begin
        nf = 6
        topo = rectangletopology(3, 2)
        obs = MajoranaSum(nf, :n, 2)
        MajoranaPropagation.pop_id!(obs)
        circ = FermionicRotation[]
        thetas = Float64[]
        for pair in topo
            push!(circ, FermionicRotation(:hop, pair))
            push!(thetas, 0.1)
        end
        for pair in topo
            push!(circ, FermionicRotation(:nn, pair))
            push!(thetas, 0.1)
        end

        obs_dict = deepcopy(obs)
        obs_vec = VectorMajoranaSum(deepcopy(obs))
        for _ in 1:2
            propagate!(circ, obs_dict, thetas; min_abs_coeff=1e-10)
            propagate!(circ, obs_vec, thetas; min_abs_coeff=1e-10)
        end

        d_dict = Dict(k => v for (k, v) in obs_dict)
        d_vec = Dict(k => v for (k, v) in obs_vec)
        @test length(d_dict) == length(d_vec)
        @test all(isapprox(v, get(d_vec, k, Inf); atol=1e-10) for (k, v) in d_dict)
    end

end
