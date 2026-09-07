#include("../src/InitialStates.jl")
#using Yao

function create_hybrid_state(nf::Int, occupied_sites, activated_links)
    nq = nf - 1
    f_part = FockState(nf, occupied_sites)
    q_part  = EigenState(nq, activated_links, :x)
    return HybridEigenState(f_part, q_part)
end 

function create_yao_state(nf::Int, occupied_sites, activated_links)
    nq = nf - 1
    state_encoding = 0
    
    for f_site in occupied_sites 
        state_encoding += 2^(2 * (f_site - 1))
    end 
    for l_site in activated_links
        state_encoding += 2^(2 * l_site - 1)
    end 

    input = zeros(ComplexF64, 2^(nf + nq))
    input[state_encoding + 1] += 1
    psi = ArrayReg(input)

    Hadamard = chain((nq + nf), put((2*ii) => Yao.H) for ii = 1:nq)
    Yao.apply!(psi, Hadamard)
    
    return psi 
end 