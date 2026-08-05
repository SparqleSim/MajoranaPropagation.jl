using MajoranaPropagation
using PauliPropagation
using Yao
using MajoranaPropagation.FermionToQubitMappings
using Test

@testset "MP vs. Jordan-Wigner PP vs. Jordan-Wigner statevector" begin
    # spinless tests 
    @testset "spinless" begin
        @testset "free fermions" begin
            n_fermions = 8
            is_spinful = false
            h = 0.2
            topo = bricklayertopology(n_fermions)

            obs = (:n, 4)
            msum = MajoranaSum(n_fermions, obs[1], obs[2])
            msum_pp = JordanWigner(msum)
            yao_obs = majoranapropagation2yao(msum)

            #initial state
            occupied_sites = 1:2:n_fermions
            fock_state = FockState(n_fermions, occupied_sites)

            yao_psi = zero_state(n_fermions)
            state_prep = chain(n_fermions, put(site => Yao.X) for site in occupied_sites)
            Yao.apply!(yao_psi, state_prep)

            #build circuits 
            circ::Vector{FermionicRotation} = []
            thetas = []
            for (i, j) in topo
                push!(circ, FermionicRotation(:hop, [i, j]))
                push!(thetas, h)
            end

            pp_circ, pp_thetas = JordanWigner(n_fermions, is_spinful, circ, thetas)
            yao_circ = majoranapropagation2yao(n_fermions, is_spinful, circ, thetas)

            mp_res = overlapwithfock(msum, fock_state)
            pp_res = overlapwithcomputational(msum_pp, occupied_sites)
            yao_res = Yao.expect(yao_obs, yao_psi)
            @test abs(mp_res - yao_res) < 1.e-12
            @test abs(pp_res - yao_res) < 1.e-12

            n_iters = 200
            for _ in 1:n_iters
                propagate!(circ, msum, thetas; min_abs_coeff=-1.)
                propagate!(pp_circ, msum_pp, pp_thetas; min_abs_coeff=-1.)
                Yao.apply!(yao_psi, yao_circ)

                mp_res = overlapwithfock(msum, fock_state)
                pp_res = overlapwithcomputational(msum_pp, occupied_sites)
                yao_res = Yao.expect(yao_obs, yao_psi)
                @test abs(mp_res - yao_res) < 1.e-12
                @test abs(pp_res - yao_res) < 1.e-12
            end
        end

        @testset "interacting fermions" begin
            n_fermions = 8
            is_spinful = false
            U = 0.5
            h = 0.2
            topo = bricklayertopology(n_fermions)

            obs = (:nn, [4, 6])
            msum = MajoranaSum(n_fermions, obs[1], obs[2])
            msum_pp = JordanWigner(msum)

            yao_obs = majoranapropagation2yao(msum)

            #initial state
            occupied_sites = 2:2:n_fermions
            fock_state = FockState(n_fermions, occupied_sites)

            yao_psi = zero_state(n_fermions)
            state_prep = chain(n_fermions, put(site => Yao.X) for site in occupied_sites)
            Yao.apply!(yao_psi, state_prep)

            #build circuits 
            circ::Vector{FermionicRotation} = []
            thetas = []
            for (i, j) in topo
                push!(circ, FermionicRotation(:hop, [i, j]))
                push!(thetas, h)
            end

            for (i, j) in topo
                push!(circ, FermionicRotation(:nn, [i, j]))
                push!(thetas, U)
            end

            pp_circ, pp_thetas = JordanWigner(n_fermions, is_spinful, circ, thetas)
            yao_circ = majoranapropagation2yao(n_fermions, is_spinful, circ, thetas)

            n_iters = 5
            for _ in 1:n_iters
                propagate!(circ, msum, thetas; min_abs_coeff=1.e-14)
                propagate!(pp_circ, msum_pp, pp_thetas; min_abs_coeff=1.e-14)
                Yao.apply!(yao_psi, yao_circ)

                mp_res = overlapwithfock(msum, fock_state)
                pp_res = overlapwithcomputational(msum_pp, occupied_sites)
                yao_res = Yao.expect(yao_obs, yao_psi)

                @test abs(mp_res - yao_res) < 1.e-12
                @test abs(pp_res - yao_res) < 1.e-12
            end
        end
    end

    # spinful tests: site i maps to qubits 2i-1 (up) and 2i (down)
    @testset "spinful" begin
        @testset "free fermions" begin
            n_sites = 4
            is_spinful = true
            nq = 2 * n_sites
            h = 0.2
            topo = bricklayertopology(n_sites)

            obs_site = 3
            msum = MajoranaSum(n_sites, :nup, obs_site)
            msum_pp = JordanWigner(msum)
            yao_obs = majoranapropagation2yao(msum)

            #initial state
            up_occupied = 1:2:n_sites
            dn_occupied = 2:2:n_sites
            fock_state = FockState(n_sites, collect(up_occupied), collect(dn_occupied))
            occupied_qubits = sort(vcat([2s - 1 for s in up_occupied], [2s for s in dn_occupied]))

            yao_psi = zero_state(nq)
            state_prep = chain(nq, put(q => Yao.X) for q in occupied_qubits)
            Yao.apply!(yao_psi, state_prep)

            #build circuits
            circ::Vector{FermionicRotation} = []
            thetas = []
            for (i, j) in topo
                push!(circ, FermionicRotation(:hopup, [i, j]))
                push!(thetas, h)
                push!(circ, FermionicRotation(:hopdn, [i, j]))
                push!(thetas, h)
            end

            pp_circ, pp_thetas = JordanWigner(n_sites, is_spinful, circ, thetas)
            yao_circ = majoranapropagation2yao(n_sites, is_spinful, circ, thetas)

            mp_res = overlapwithfock(msum, fock_state)
            pp_res = overlapwithcomputational(msum_pp, occupied_qubits)
            yao_res = Yao.expect(yao_obs, yao_psi)
            @test abs(mp_res - yao_res) < 1.e-12
            @test abs(pp_res - yao_res) < 1.e-12

            n_iters = 200
            for _ in 1:n_iters
                propagate!(circ, msum, thetas; min_abs_coeff=-1.)
                propagate!(pp_circ, msum_pp, pp_thetas; min_abs_coeff=-1.)
                Yao.apply!(yao_psi, yao_circ)

                mp_res = overlapwithfock(msum, fock_state)
                pp_res = overlapwithcomputational(msum_pp, occupied_qubits)
                yao_res = Yao.expect(yao_obs, yao_psi)
                @test abs(mp_res - yao_res) < 1.e-12
                @test abs(pp_res - yao_res) < 1.e-12
            end
        end

        @testset "Hubbard" begin
            n_sites = 4
            is_spinful = true
            nq = 2 * n_sites
            U = 0.5
            h = 0.2
            topo = bricklayertopology(n_sites)

            obs_site = 3
            msum = MajoranaSum(n_sites, :nupndn, obs_site)
            msum_pp = JordanWigner(msum)

            yao_obs = majoranapropagation2yao(msum)

            #initial state
            up_occupied = 1:2:n_sites
            dn_occupied = 2:2:n_sites
            fock_state = FockState(n_sites, collect(up_occupied), collect(dn_occupied))
            occupied_qubits = sort(vcat([2s - 1 for s in up_occupied], [2s for s in dn_occupied]))

            yao_psi = zero_state(nq)
            state_prep = chain(nq, put(q => Yao.X) for q in occupied_qubits)
            Yao.apply!(yao_psi, state_prep)

            #build circuits
            circ::Vector{FermionicRotation} = []
            thetas = []
            for (i, j) in topo
                push!(circ, FermionicRotation(:hopup, [i, j]))
                push!(thetas, h)
                push!(circ, FermionicRotation(:hopdn, [i, j]))
                push!(thetas, h)
            end
            for i = 1:n_sites
                push!(circ, FermionicRotation(:nupndn, i))
                push!(thetas, U)
            end

            pp_circ, pp_thetas = JordanWigner(n_sites, is_spinful, circ, thetas)
            yao_circ = majoranapropagation2yao(n_sites, is_spinful, circ, thetas)

            n_iters = 5
            for _ in 1:n_iters
                propagate!(circ, msum, thetas; min_abs_coeff=1.e-14)
                propagate!(pp_circ, msum_pp, pp_thetas; min_abs_coeff=1.e-14)
                Yao.apply!(yao_psi, yao_circ)

                mp_res = overlapwithfock(msum, fock_state)
                pp_res = overlapwithcomputational(msum_pp, occupied_qubits)
                yao_res = Yao.expect(yao_obs, yao_psi)

                @test abs(mp_res - yao_res) < 1.e-12
                @test abs(pp_res - yao_res) < 1.e-12
            end
        end
    end
end