using Yao

include("yao_helpers/fermionicgates_to_yao.jl")
@testset "Majorana algebra" begin
    # spinless tests 
    @testset "spinless" begin
        @testset "free fermions" begin
            n_fermions = 12
            h = 0.2
            topo = bricklayertopology(n_fermions)

            obs = (:n, 4)
            msum = MajoranaSum(n_fermions, obs[1], obs[2])
            yao_obs = kron(n_fermions, obs[2] => Z)

            #initial state
            occupied_sites = 1:2:n_fermions
            fock_state = FockState(n_fermions, occupied_sites)

            yao_psi = zero_state(n_fermions)
            state_prep = chain(n_fermions, put(site => Yao.X) for site in occupied_sites)
            Yao.apply!(yao_psi, state_prep)

            #build circuits 
            circ::Vector{FermionicGate} = []
            thetas = []
            for (i, j) in topo
                push!(circ, FermionicGate(:hop, [i, j]))
                push!(thetas, h)
            end

            yao_circ = circ_to_yao(n_fermions, circ, thetas)

            n_iters = 200
            for _ in 1:n_iters
                propagate!(circ, msum, thetas; min_abs_coeff=-1.)
                Yao.apply!(yao_psi, yao_circ)

                mp_res = overlapwithfock(msum, fock_state)
                yao_res = (1. - Yao.expect(yao_obs, yao_psi)) / 2. # convert from Z expectation to occupation number
                #@show mp_res, yao_res
                @test abs(mp_res - yao_res) < 1.e-12
            end
        end

        @testset "interacting fermions" begin
            n_fermions = 10
            U = 0.5
            h = 0.2
            topo = bricklayertopology(n_fermions)

            obs = (:nn, [4, 6])
            msum = MajoranaSum(n_fermions, obs[1], obs[2])
            
            yao_Z_s1 = kron(n_fermions, obs[2][1] => Z)
            yao_Z_s2 = kron(n_fermions, obs[2][2] => Z)
            yao_ZZ = kron(n_fermions, site => Z for site in obs[2])

            #initial state
            occupied_sites = 2:2:n_fermions
            fock_state = FockState(n_fermions, occupied_sites)

            yao_psi = zero_state(n_fermions)
            state_prep = chain(n_fermions, put(site => Yao.X) for site in occupied_sites)
            Yao.apply!(yao_psi, state_prep)

            #build circuits 
            circ::Vector{FermionicGate} = []
            thetas = []
            for (i, j) in topo
                push!(circ, FermionicGate(:hop, [i, j]))
                push!(thetas, h)
            end

            for (i, j) in topo
                push!(circ, FermionicGate(:nn, [i, j]))
                push!(thetas, U)
            end

            yao_circ = circ_to_yao(n_fermions, circ, thetas)

            n_iters = 5
            for _ in 1:n_iters
                propagate!(circ, msum, thetas; min_abs_coeff=1.e-14)
                Yao.apply!(yao_psi, yao_circ)

                mp_res = overlapwithfock(msum, fock_state)
                yao_res = (1. - Yao.expect(yao_Z_s1, yao_psi) - Yao.expect(yao_Z_s2, yao_psi) + Yao.expect(yao_ZZ, yao_psi)) / 4. # convert from Z expectations to occupation number
                #@show mp_res, yao_res
                @test abs(mp_res - yao_res) < 1.e-12
            end
        end
    end
end