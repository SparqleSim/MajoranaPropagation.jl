using MajoranaPropagation
using PauliPropagation

using Test
using Random
Random.seed!(42)

function create_omega_L(nfermions)
    omega_L = zeros(Int, 2 * nfermions, 2 * nfermions)
    for i = 1:2*nfermions
        for j = 1:i-1
            omega_L[i, j] = 1
        end
    end
    return omega_L
end

function create_omega(nfermions)
    omega = create_omega_L(nfermions)
    omega .+= transpose(omega)
    return omega
end

function compute_commutation(ms1::Vector{Int}, ms2::Vector{Int}, omega)
    return (-1)^mod(ms1' * omega * ms2, 2)
end

function compute_prefactor(ms1::Vector{Int}, ms2::Vector{Int}, omega, omega_L)
    g_val = ms1' * omega_L * ms2 + (ms1' * omega_L * ms1) * (ms2' * omega_L * ms2) + (ms1' * omega * ms2) * ((ms1' * omega_L * ms1) + (ms2' * omega_L * ms2) + 1)
    return (-1)^mod(g_val, 2) * (1im)^mod(ms1' * omega * ms2, 2)
end

function gamma_to_dense(ms::MajoranaString)
    n = ms.nfermions
    dense_vec = zeros(Int, 2 * n)
    for i in 1:2*n
        if (ms.gammas >> (i-1)) & 1 == 1
            dense_vec[i] = 1
        end
    end
    return dense_vec
end


@testset "Commutation relations" begin
    @testset "single site" begin
        nf = 1
        omega_L_dense = create_omega_L(nf)
        omega_dense = create_omega(nf)
        gamma = MajoranaString(nf, [1])
        gamma_prime = MajoranaString(nf, [2])
        density = MajoranaString(nf, [1, 2])

        gamma_dense = gamma_to_dense(gamma)
        gamma_prime_dense = gamma_to_dense(gamma_prime)
        density_dense = gamma_to_dense(density)

        @test omega_mult(gamma, gamma_prime) == mod(gamma_dense' * omega_dense * gamma_prime_dense, 2)
        @test omega_mult(gamma, density) == mod(density_dense' * omega_dense * gamma_dense, 2)
    end

    @testset "multi-sites" begin
        nf = 20
        omega_L_dense = create_omega_L(nf)
        omega_dense = create_omega(nf)

        gamma1 = MajoranaString(nf, [1])
        gamma_prime1 = MajoranaString(nf, [2])
        density1 = MajoranaString(nf, [1, 2])
        density2 = MajoranaString(nf, [3, 4])
        gamma2 = MajoranaString(nf, [3])
        gamma_prime2 = MajoranaString(nf, [4])
        gamma3 = MajoranaString(nf, [5])
        gamma_prime3 = MajoranaString(nf, [6])

        gamma1_dense = gamma_to_dense(gamma1)
        gamma_prime1_dense = gamma_to_dense(gamma_prime1)
        density1_dense = gamma_to_dense(density1)
        density2_dense = gamma_to_dense(density2)
        gamma2_dense = gamma_to_dense(gamma2)
        gamma_prime2_dense = gamma_to_dense(gamma_prime2)
        gamma3_dense = gamma_to_dense(gamma3)
        gamma_prime3_dense = gamma_to_dense(gamma_prime3)

        #check products
        pref, g1p3 = ms_mult(gamma1, gamma_prime3)
        pref, p1g2 = ms_mult(gamma_prime1, gamma2)
        pref, g1p2 = ms_mult(gamma1, gamma_prime2)

        g1p3_dense = gamma_to_dense(g1p3)
        p1g2_dense = gamma_to_dense(p1g2)
        g1p2_dense = gamma_to_dense(g1p2)

        @test pref == compute_prefactor(gamma1_dense, gamma_prime3_dense, omega_dense, omega_L_dense)
        @test pref == compute_prefactor(gamma_prime1_dense, gamma2_dense, omega_dense, omega_L_dense)

        #check single site test in multi sites
        @test omega_mult(gamma1, gamma_prime1) == mod(gamma1_dense' * omega_dense * gamma_prime1_dense, 2)
        @test omega_mult(density1, gamma_prime1) == mod(density1_dense' * omega_dense * gamma_prime1_dense, 2)

        @test omega_mult(density1, density2) == mod(density1_dense' * omega_dense * density2_dense, 2)
        @test omega_mult(gamma1, density2) == mod(gamma1_dense' * omega_dense * density2_dense, 2)
        @test omega_mult(gamma_prime1, density2) == mod(gamma_prime1_dense' * omega_dense * density2_dense, 2)

        @test omega_mult(g1p3, p1g2) == mod(g1p3_dense' * omega_dense * p1g2_dense, 2)
        @test omega_mult(g1p3, gamma2) == mod(g1p3_dense' * omega_dense * gamma2_dense, 2)
        @test omega_mult(g1p3, gamma1) == mod(g1p3_dense' * omega_dense * gamma1_dense, 2)
        @test omega_mult(g1p3, g1p2) == mod(g1p3_dense' * omega_dense * g1p2_dense, 2)


        ms1 = MajoranaString(nf, [1, 4, 6, 10, 11])
        ms2 = MajoranaString(nf, [2, 15])
        ms3 = MajoranaString(nf, [3, 5, 7, 8, 9])

        ms1_dense = gamma_to_dense(ms1)
        ms2_dense = gamma_to_dense(ms2)
        ms3_dense = gamma_to_dense(ms3)

        @test omega_mult(ms1, ms2) == mod(ms1_dense' * omega_dense * ms2_dense, 2)
        @test omega_mult(ms1, ms3) == mod(ms1_dense' * omega_dense * ms3_dense, 2)

    end

    @testset "multi-sites small, arbitrary strings" begin
        nf = 10

        omega_dense = create_omega(nf)
        omega_L_dense = create_omega_L(nf)
        
        TT = getinttype(nf)
        max_val = MajoranaString(nf, [2 * nf]).gammas

        n_tests = 2000
        for _ = 1:n_tests
            ms1 = rand(TT)
            while get_weight(ms1) % 2 != 0 || ms1 > max_val
                ms1 = rand(TT)
            end
            ms2 = rand(TT)
            while get_weight(ms2) % 2 != 0 || ms2 > max_val
                ms2 = rand(TT)
            end

            ms1 = MajoranaString(nf, ms1)
            ms2 = MajoranaString(nf, ms2)

            ms1_dense = gamma_to_dense(ms1)
            ms2_dense = gamma_to_dense(ms2)

            @test omega_mult(ms1, ms2) == mod(ms1_dense' * omega_dense * ms2_dense, 2)
            @test omega_L_mult(ms1, ms2) == mod(ms1_dense' * omega_L_dense * ms2_dense, 2)
            @test omega_L_mult(ms1) == mod(ms1_dense' * omega_L_dense * ms1_dense, 2)
            @test omega_L_mult(ms2) == mod(ms2_dense' * omega_L_dense * ms2_dense, 2)

            pref, _ = ms_mult(ms1, ms2)

            @test pref == compute_prefactor(ms1_dense, ms2_dense, omega_dense, omega_L_dense)
        end

    end

    @testset "multi-sites large, arbitrary strings" begin
        nf = 100

        omega_dense = create_omega(nf)
        omega_L_dense = create_omega_L(nf)

        TT = getinttype(nf)
        max_val = MajoranaString(nf, [2 * nf]).gammas

        n_tests = 2000
        for _ = 1:n_tests
            ms1 = rand(TT)
            while get_weight(ms1) % 2 != 0 || ms1 > max_val
                ms1 = rand(TT)
            end
            ms2 = rand(TT)
            while get_weight(ms2) % 2 != 0 || ms2 > max_val
                ms2 = rand(TT)
            end

            ms1 = MajoranaString(nf, ms1)
            ms2 = MajoranaString(nf, ms2)

            ms1_dense = gamma_to_dense(ms1)
            ms2_dense = gamma_to_dense(ms2)

            @test omega_mult(ms1, ms2) == mod(ms1_dense' * omega_dense * ms2_dense, 2)
            @test omega_L_mult(ms1, ms2) == mod(ms1_dense' * omega_L_dense * ms2_dense, 2)
            @test omega_L_mult(ms1) == mod(ms1_dense' * omega_L_dense * ms1_dense, 2)
            @test omega_L_mult(ms2) == mod(ms2_dense' * omega_L_dense * ms2_dense, 2)

            pref, _ = ms_mult(ms1, ms2)

            @test pref == compute_prefactor(ms1_dense, ms2_dense, omega_dense, omega_L_dense)
        end

    end

end

@testset "Gate applications" begin

    @testset "small" begin
        nfermions = 2
        omega_dense = create_omega(nfermions)
        omega_L_dense = create_omega_L(nfermions)

        mu_g = MajoranaString(nfermions, [1, 4])
        mu_v = MajoranaString(nfermions, [1, 2])
        gate = MajoranaRotation(mu_g)

        msum = MajoranaSum(Float64, nfermions, false)
        PauliPropagation.PropagationBase.add!(msum, mu_v.gammas, 1.)
        theta = 0.5
        propagate!([gate], msum, [theta])

        mu_g_dense = gamma_to_dense(mu_g)
        mu_v_dense = gamma_to_dense(mu_v)

        pref = compute_prefactor(mu_g_dense, mu_v_dense, omega_dense, omega_L_dense)
        _, ms_prod = ms_mult(mu_g, mu_v)

        msum_analytic = MajoranaSum(Float64, nfermions, false)
        PauliPropagation.PropagationBase.add!(msum_analytic, mu_v.gammas, cos(theta))
        PauliPropagation.PropagationBase.add!(msum_analytic, ms_prod.gammas, real(sin(theta) * 1im * pref))
        @test msum == msum_analytic
    end

    @testset "large" begin
        nf = 50
        omega_dense = create_omega(nf)
        omega_L_dense = create_omega_L(nf)
        
        TT = getinttype(nf)
        max_val = MajoranaString(nf, [2 * nf]).gammas

        n_tests = 2000
        for _ = 1:n_tests
            ms_v = rand(TT)
            while get_weight(ms_v) % 2 != 0 || ms_v > max_val
                ms_v = rand(TT)
            end
            ms_gate = rand(TT)
            while get_weight(ms_gate) % 2 != 0 || ms_gate > max_val || MajoranaPropagation.commutes(ms_gate, ms_v)
                ms_gate = rand(TT)
            end

            ms_v = MajoranaString(nf, ms_v)
            ms_gate = MajoranaString(nf, ms_gate)

            ms_v_dense = gamma_to_dense(ms_v)
            ms_gate_dense = gamma_to_dense(ms_gate)

            pref = compute_prefactor(ms_gate_dense, ms_v_dense, omega_dense, omega_L_dense)
            _, ms_prod = ms_mult(ms_gate, ms_v)

            theta_angle = rand() * 2 * π

            msum = MajoranaSum(Float64, nf, false)
            PauliPropagation.PropagationBase.add!(msum, ms_v.gammas, 1.)
            propagate!([MajoranaRotation(ms_gate)], msum, [theta_angle])

            msum_analytic = MajoranaSum(Float64, nf, false)
            PauliPropagation.PropagationBase.add!(msum_analytic, ms_v.gammas, cos(theta_angle))
            sin_pref = real(sin(theta_angle) * 1im * pref)
            PauliPropagation.PropagationBase.add!(msum_analytic, ms_prod.gammas, sin_pref)

            @test msum == msum_analytic
        end
    end
end