using ITensors, ITensorMPS

function mpo_circ(circ, sites, thetas)
    res::Vector{ITensor} = []
    for (j, gate) in enumerate(circ)
        if gate[1] == "n"
            push!(res, exp(-im*thetas[j]*op("N",sites[gate[2]])))
        elseif gate[1] == "h"
            push!(res, exp(-im*thetas[j]*(op("Cdag",sites[gate[2]])*op("C",sites[gate[3]])+op("Cdag",sites[gate[3]])*op("C",sites[gate[2]]))))
        elseif gate[1] == "nund"
            push!(res, exp(-im*thetas[j]*op("N",sites[gate[2]])*op("N",sites[gate[3]])))
        else 
            error("Gate $(gate) not recognized")
        end
    end 
    return res
end

function build_psi(initial_state, circ_mpo, sites, N, χ)
    initial_state_explicity = []
    for j=1:N
        push!(initial_state_explicity, j in initial_state ? "1" : "0")
    end
    psi0 =  MPS(sites, initial_state_explicity)
    psi = ITensorMPS.apply(circ_mpo, psi0; maxdim=χ)
    return psi
end

function mpo_obs(identifier::String, acting_on, sites, N)
    mpo_obs = MPO(sites)
    for j in 1:N
        if j in acting_on
            mpo_obs[j] = op(identifier, sites[j])
        else
            mpo_obs[j] = op("Id", sites[j])
        end
    end
    return mpo_obs
end 


function mpo_obs(identifiers, acting_on, sites, N)
    mpo_obs = MPO(sites)
    for j in 1:N
        mpo_obs[j] = op("Id", sites[j])
    end

    for (j, identifier) in enumerate(identifiers)
        mpo_obs[acting_on[j]] = op(identifier, sites[acting_on[j]])
    end
    return mpo_obs
end 

function mpo_hopping(acting_on, sites, N)
    #error("Not implemented")
    o1 = MPO(sites)
    o2 = MPO(sites)
    for j in 1:N
        if j==acting_on[1]
            o1[j] = op("Cdag", sites[j])
            o2[j] = op("C", sites[j])
        elseif  j==acting_on[2]
            o1[j] = op("C", sites[j])
            o2[j] = op("Cdag", sites[j])
        else 
            o1[j] = op("Id", sites[j])
            o2[j] = op("Id", sites[j])
        end 
    end
    return o1 + o2
end

function run_circ_mps(N, circ, initial_state, thetas, list_of_obs::String; χ=150)
    @assert list_of_obs == "N"
    sites = siteinds("Fermion",N)
    circ_mpo = mpo_circ(circ, sites, thetas)
    psi = build_psi(initial_state, circ_mpo, sites, N, χ)
    return expect(psi,"N")
end 
    
function run_circ_mps(N, circ, initial_state, thetas, list_of_obs; χ=150)
    sites = siteinds("Fermion",N)
    circ_mpo = mpo_circ(circ, sites, thetas)
    psi = build_psi(initial_state, circ_mpo, sites, N, χ)

    obs_res = zeros(length(list_of_obs))
    for (j, obs) in enumerate(list_of_obs)
        if obs[1] == "cdag_c"
            mpo_obs_loc = mpo_hopping(obs[2], sites, N)
        else
            mpo_obs_loc = mpo_obs(obs[1], obs[2], sites, N)
        end
        obs_res[j] = real(inner(psi', mpo_obs_loc, psi))
    end
    return obs_res
end 