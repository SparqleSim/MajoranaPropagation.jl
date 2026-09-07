#using Yao
#include("../src/HybridSum.jl")

function create_yao_obs(n_fermions::Int, n_i::Int=0, X_is::Array{Int}=Int[])
    n_tot = 2 * n_fermions - 1
    X_is *= 2
    n_i = 2 * n_i - 1

    obs = kron(n_tot, i=>X for i in X_is)
    if n_i > 0
        push!(X_is, n_i)
        obs = 0.5 * obs - 0.5 * kron(n_tot, i=>(0 == mod(i,2) ? X : Z) for i in X_is)
    end 
    return obs
end 

function create_hybrid_obs(n_fermions::Int, n_i::Int=0, X_is::Array{Int}=Int[])
    n_qubits = n_fermions - 1
    symbs = [:X for i in X_is]
    msum = MajoranaSum(n_fermions, :n, n_i)
    if n_i == 0
        msum = MajoranaSum(n_fermions)
    end
    pstr = PauliString(n_qubits, symbs, X_is)
    return HybridSum(msum, pstr)
end 