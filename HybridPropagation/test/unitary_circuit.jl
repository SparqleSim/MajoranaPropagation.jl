#using Yao
#include("../src/HybridGates.jl")

#---Yao Circuit-------------------------------------------------
function single_ZZZ_term!(circ, finds, qind, theta::Float64)
    push!(circ, cnot(finds[1], qind))
    push!(circ, cnot(finds[2], qind))
    push!(circ, Yao.put(qind=>Yao.Rz(theta)))
    push!(circ, cnot(finds[2], qind))
    push!(circ, cnot(finds[1], qind))
    return circ 
end 

function single_XXZ_term!(circ, f_ind, q_ind::Int, theta::Float64)
    for i in 1:2
        push!(circ, Yao.put(f_ind[i] => Yao.H))
    end 
    single_ZZZ_term!(circ, f_ind, q_ind, theta)
    for i in 1:2
        push!(circ, Yao.put(f_ind[i] => Yao.H))
    end 
    return circ
end 

function single_YYZ_term!(circ, f_ind, q_ind::Int, theta::Float64)
    for i in 1:2
        push!(circ, Yao.put(f_ind[i] => Yao.shift(-π/2)))
        push!(circ, Yao.put(f_ind[i] => Yao.H))
    end 
    single_ZZZ_term!(circ, f_ind, q_ind, theta)
    for i in 1:2
        push!(circ, Yao.put(f_ind[i] => Yao.H))
        push!(circ, Yao.put(f_ind[i] => Yao.shift(π/2)))
    end 
    return circ
end 

function single_X_term!(circ, qind::Int, theta::Float64)
    push!(circ, Yao.put(qind => Yao.Rx(theta)))
    return circ
end

function create_circuit(n_fermionic_sites::Int, del_t::Float64; t::Float64=1.0, h::Float64=1.0)
    n_links = n_fermionic_sites - 1
    n_tot = n_fermionic_sites + n_links

    circuit = chain(n_tot)
    link_site_ind = [(2 * ind) for ind in 1:n_links]
    for link in link_site_ind
        single_XXZ_term!(circuit, ((link-1), (link+1)), link, (-t*del_t)) #2 * 0.5
        single_YYZ_term!(circuit, ((link-1), (link+1)), link, (-t*del_t)) #2 * 0.5
    end

    for link in link_site_ind
        single_X_term!(circuit, link, (-2*h*del_t))
    end 

    return circuit
end


#---Hybrid Unitary----------------------------------------
function create_unitary(n_fermionic_sites::Int, del_t::Float64; t::Float64=1.0, h::Float64=1.0)
    n_links = n_fermionic_sites - 1
    circuit = Gate[]
    for link in 1:n_links
        push!(circuit, HybridGate(n_fermionic_sites, :hop, [link, (link+1)], n_links, [:Z], [link], false))
    end 

    for link in 1:n_links
        push!(circuit, HybridGate(n_fermionic_sites, n_links, [:X], [link]))
    end 

    thetas = ones(Float64, (2 * n_links))
    thetas[1:n_links] *= -t * del_t
    thetas[(n_links+1):end] *= -h * del_t

    return circuit, thetas
end 