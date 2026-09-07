using Test
using ProgressMeter
using Yao

include("../src/HybridPropagation.jl")

include("unitary_circuit.jl")
include("state_preparation.jl")
include("obs_preparation.jl")

eps = 1.e-12    #Defines comparison threshold

@testset "Hybrid Propagation" begin
    
    @testset "Link Measurement" begin 
        #System Definition
        n_fermions = 10
        h = 1.2
        t = 0.5

        #Simulation Definition
        del_t = 0.1
        nTSteps = 20

        #Observable = X_middle
        X_middle = Int[Int(n_fermions/2)]
        n_i = 0

        yao_obs = create_yao_obs(n_fermions, n_i, X_middle)
        hybrid_obs = create_hybrid_obs(n_fermions, n_i, X_middle)

        #Initial State
        initial_occupied_fermionic_sites = [i for i in 1:2:n_fermions]
        initial_active_links = []

        yao_state = create_yao_state(n_fermions, initial_occupied_fermionic_sites, initial_active_links)
        hybrid_state = create_hybrid_state(n_fermions, initial_occupied_fermionic_sites, initial_active_links)

        #Circuit
        yao_circ = create_circuit(n_fermions, del_t)
        hybrid_circ, thetas = create_unitary(n_fermions, del_t)

        f_filter, q_filter = create_filters(hybrid_obs)
        @showprogress for t in 1:nTSteps
            propagate!(hybrid_circ, hybrid_obs, thetas; fermions_filter=f_filter, qubits_filter = q_filter, min_abs_coeff=-1.)
            hybrid_expval = overlapwithstate(hybrid_obs, hybrid_state)

            Yao.apply!(yao_state, yao_circ)
            yao_expval = Yao.expect(yao_obs, yao_state)

            #@show hybrid_expval, yao_expval
            @test abs(yao_expval- hybrid_expval) < eps
        end 
    end 
    
    @testset "Fermion Number Measurement" begin 
        #System Definition
        n_fermions = 10
        h = 0.7
        t = 1.3

        #Simulation Definition
        del_t = 0.1
        nTSteps = 20

        #Observable = n_5
        n_i = 5
        X_i = Int[]

        yao_obs = create_yao_obs(n_fermions, n_i, X_i)
        hybrid_obs = create_hybrid_obs(n_fermions, n_i, X_i)

        #Initial State
        initial_occupied_fermionic_sites = [i for i in 1:2:n_fermions]
        initial_active_links = []

        yao_state = create_yao_state(n_fermions, initial_occupied_fermionic_sites, initial_active_links)
        hybrid_state = create_hybrid_state(n_fermions, initial_occupied_fermionic_sites, initial_active_links)

        #Circuit
        yao_circ = create_circuit(n_fermions, del_t)
        hybrid_circ, thetas = create_unitary(n_fermions, del_t)

        f_filter, q_filter = create_filters(hybrid_obs)
        @showprogress for t in 1:nTSteps
            propagate!(hybrid_circ, hybrid_obs, thetas; fermions_filter=f_filter, qubits_filter = q_filter, min_abs_coeff=-1.)
            hybrid_expval = overlapwithstate(hybrid_obs, hybrid_state)

            Yao.apply!(yao_state, yao_circ)
            yao_expval = Yao.expect(yao_obs, yao_state)

            #@show hybrid_expval, yao_expval
            @test abs(yao_expval- hybrid_expval) < eps
        end 
    end
    
end 