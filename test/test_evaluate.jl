using Test
using Random
using ITensors, ITensorMPS
Random.seed!(42)

#@testset "Evaluate" begin
    
#end

n_fermions = 4 

m = [1, 4]
n = [2, 3]

ms1 = MajoranaString(n_fermions, [2*1-1, 2*2])

res = fockevaluate(ms1, m, n)
@show res 

#mps stuff 
sites = siteinds("Fermion", n_fermions)
m_expl = []
n_expl = []
for j=1:n_fermions
    push!(m_expl, j in m ? "1" : "0")
    push!(n_expl, j in n ? "1" : "0")
end
n_psi =  MPS(sites, n_expl)
m_psi =  MPS(sites, m_expl)

ms_op = MPO(sites, "Id")

for i=1:n_fermions
    ops = ["Id" for _ in 1:n_fermions]
    gamma_2i = 1im * (op("Cdag", sites[i]) - op("C", sites[i]))
    ops_tensors = [op("Id", sites[j]) for j in 1:n_fermions]
    ops_tensors[i] = gamma_2i
    majorana_mpo = MPO(ops_tensors, sites)
    #=@show typeof(op("Cdag", sites[i]))
    @show typeof(op(sites, "Cdag", i))
    @show MPO(op(sites, "Cdag", i), sites)
    if ((ms1.gammas >> (2*i-2)) & 1) == 1
        
        #global ms_op *= MPO((op("Cdag", sites[i]) + op("C", sites[i])), sites)
    end=#
end
@show ms_op
