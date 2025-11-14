""" 
    general_trotter_layer(groups::Vector{Vector{Tuple{FermionicGate,Float64}}}, group_indices::Vector{Int}, invert_groups::Vector{Bool}, group_prefactors::Vector{Float64})
General implementation of a single Trotter layer.
Given the terms in the Hamiltonian are separated in `groups`, implement a Trotter step 
exp(-i groups[j_1] * group_prefactors[j_1]) ... exp(groups[j_N] * group_prefactors[j_N])
where the order of appearence j_1, ..., j_N = group_indices (a single group can appear multiple times in higher order trotterizations),
and different Trotterization schemes can be specified by changing `group_prefactors`.
The boolean `invert_groups[j]` prescribes the current group is to be implemented in reverse oder or not.
"""
function general_trotter_layer(groups::Vector{Vector{Tuple{FermionicGate,Float64}}}, group_indices::Vector{Int}, invert_groups::Vector{Bool}, group_prefactors::Vector{Float64})
    @assert length(group_indices) == length(invert_groups)
    @assert length(group_indices) == length(group_prefactors)

    circuit::Vector{FermionicGate} = []
    thetas::Vector{Float64} = []

    for j = 1:length(group_indices)
        group_index = group_indices[j]
        group_to_pass = invert_groups[j] == false ? groups[group_index] : reverse(groups[group_index])
        for (gate, theta) in group_to_pass
            push!(circuit, gate)
            push!(thetas, theta * group_prefactors[j])
        end
    end

    return circuit, thetas
end


""" 
    trotter_layer(groups::Vector{Vector{Tuple{FermionicGate,Float64}}}, order::Symbol; rotation_prefactor=2., kwargs...)
Construct a trotter layer for a specified `order` from groups of fermionic gates and their associated angles.
Currently supported orders are `:first` and `:second`.

`:first` implements a first-order Trotter layer by applying each group of gates sequentially.
`:second` implements a second-order Trotter layer by applying half of the gates from each group in the forward order, then applying the last group fully, and finally applying the first groups in reverse order with half angles.

`rotation_prefactor` is set to 2. by default as the Majorana rotations implement exp(-i * theta/2 * mstring), whereas fermionic gates are typically defined as exp(-i * theta * operator). This prefactor adjusts the angles accordingly.
"""
function trotter_layer(groups::Vector{Vector{Tuple{FermionicGate,Float64}}}, order::Symbol; rotation_prefactor=2., kwargs...)
    return trotter_layer(groups, Val(order); rotation_prefactor=rotation_prefactor, kwargs...)
end

""" 
    first order trotter layer
"""
function trotter_layer(groups::Vector{Vector{Tuple{FermionicGate,Float64}}}, ::Val{:first}; rotation_prefactor=2., kwargs...)
    group_indices::Vector{Int} = collect(1:length(groups))
    invert_groups::Vector{Bool} = fill(false, length(groups))
    group_prefactors::Vector{Float64} = fill(rotation_prefactor, length(groups))
    return general_trotter_layer(groups, group_indices, invert_groups, group_prefactors)
end

""" 
    second order trotter layer
"""
function trotter_layer(groups::Vector{Vector{Tuple{FermionicGate,Float64}}}, ::Val{:second}; rotation_prefactor=2., kwargs...)
    group_indices::Vector{Int} = vcat(collect(1:length(groups)), reverse(collect(1:length(groups)-1)))
    invert_groups::Vector{Bool} = vcat(fill(false, length(groups)), fill(true, length(groups) - 1))
    group_prefactors::Vector{Float64} = vcat(fill(rotation_prefactor / 2., length(groups) - 1),
        [rotation_prefactor],
        fill(rotation_prefactor / 2., length(groups) - 1))
    return general_trotter_layer(groups, group_indices, invert_groups, group_prefactors)
end

function trotter_layer(groups::Vector{Vector{Tuple{FermionicGate,Float64}}}, ::Val{symb}; kwargs...) where {symb}
    error("Trotter order $symb not recognized.")
end

function hubbard_circ_fermionic_sites_single_layer(topology, N_spinful_sites::Int, t::Float64, U::Float64, dt::Float64; return_mps_instructions=false, return_separated=false)
    mps_instructions = []
    mps_thetas = []

    #down part
    circ_down_hopping::Vector{MajoranaRotation} = []
    thetas_down_hopping::Vector{Float64} = []
    for (i, j) in topology
        @assert i < j
        #@show i, j
        ms_hop_term_i_jprime = MajoranaString(2 * N_spinful_sites, [4 * i - 1, 4 * j])
        ms_hop_term_iprime_j = MajoranaString(2 * N_spinful_sites, [4 * i, 4 * j - 1])

        push!(circ_down_hopping, MajoranaRotation(ms_hop_term_i_jprime))
        push!(thetas_down_hopping, -t * dt)
        push!(circ_down_hopping, MajoranaRotation(ms_hop_term_iprime_j))
        push!(thetas_down_hopping, +t * dt)

        if return_mps_instructions
            push!(mps_instructions, ("h", 2 * i, 2 * j))
            push!(mps_thetas, -dt * t)
        end
    end

    #up part
    circ_up_hopping::Vector{MajoranaRotation} = []
    thetas_up_hopping::Vector{Float64} = []
    for (i, j) in topology
        @assert i < j
        ms_hop_term_i_jprime = MajoranaString(2 * N_spinful_sites, [4 * i - 3, 4 * j - 2])
        ms_hop_term_iprime_j = MajoranaString(2 * N_spinful_sites, [4 * i - 2, 4 * j - 3])

        push!(circ_up_hopping, MajoranaRotation(ms_hop_term_i_jprime))
        push!(thetas_up_hopping, -t * dt)
        push!(circ_up_hopping, MajoranaRotation(ms_hop_term_iprime_j))
        push!(thetas_up_hopping, +t * dt)

        if return_mps_instructions
            push!(mps_instructions, ("h", 2 * i - 1, 2 * j - 1))
            push!(mps_thetas, -dt * t)
        end
    end

    #repulsion term
    circ_repulsion::Vector{MajoranaRotation} = []
    thetas_repulsion::Vector{Float64} = []
    for k = 1:N_spinful_sites
        ms_num_term_up = MajoranaString(2 * N_spinful_sites, [4 * k - 3, 4 * k - 2])
        ms_num_term_down = MajoranaString(2 * N_spinful_sites, [4 * k - 1, 4 * k])
        #@show ms_num_term_up, ms_num_term_down
        up_down_pref, ms_num_term_up_down = ms_mult(ms_num_term_up, ms_num_term_down)
        #@show up_down_pref, ms_num_term_up_down
        @assert imag(up_down_pref) < 1.e-12
        up_down_pref = real(up_down_pref)

        push!(circ_repulsion, MajoranaRotation(ms_num_term_up))
        push!(thetas_repulsion, U * dt / 2)
        push!(circ_repulsion, MajoranaRotation(ms_num_term_down))
        push!(thetas_repulsion, U * dt / 2)
        push!(circ_repulsion, MajoranaRotation(ms_num_term_up_down))
        push!(thetas_repulsion, up_down_pref * U * dt / 2)

        if return_mps_instructions
            push!(mps_instructions, ("nund", 2 * k - 1, 2 * k))
            push!(mps_thetas, dt * U)
        end

    end

    if return_separated
        circs = [circ_down_hopping, circ_up_hopping, circ_repulsion]
        thetas = [thetas_down_hopping, thetas_up_hopping, thetas_repulsion]
    else
        circs = vcat(circ_down_hopping, circ_up_hopping, circ_repulsion)
        thetas = vcat(thetas_down_hopping, thetas_up_hopping, thetas_repulsion)
    end

    if return_mps_instructions
        return circs, thetas, mps_instructions, mps_thetas
    end
    return circs, thetas
end

function fermionic_hubbard_circ_fermionic_sites_single_layer(topology, N_spinful_sites::Int, t::Float64, U::Float64, dt::Float64; return_mps_instructions=false, return_separated=false)
    mps_instructions = []
    mps_thetas = []

    #down part
    circ_down_hopping::Vector{FermionicGate} = []
    thetas_down_hopping::Vector{Float64} = []
    for (i, j) in topology
        @assert i < j
        #@show i, j
        hopping_ms::Vector{MajoranaString} = []
        hopping_pref::Vector{Float64} = []
        ms_hop_term_i_jprime = MajoranaString(2 * N_spinful_sites, [4 * i - 1, 4 * j])
        ms_hop_term_iprime_j = MajoranaString(2 * N_spinful_sites, [4 * i, 4 * j - 1])
        push!(hopping_ms, ms_hop_term_i_jprime)
        push!(hopping_ms, ms_hop_term_iprime_j)
        push!(hopping_pref, +1)
        push!(hopping_pref, -1)

        push!(circ_down_hopping, FermionicGate(hopping_ms, hopping_pref))
        push!(thetas_down_hopping, -t * dt)

        if return_mps_instructions
            push!(mps_instructions, ("h", 2 * i, 2 * j))
            push!(mps_thetas, -dt * t)
        end
    end

    #up part
    circ_up_hopping::Vector{FermionicGate} = []
    thetas_up_hopping::Vector{Float64} = []
    for (i, j) in topology
        @assert i < j
        hopping_ms::Vector{MajoranaString} = []
        hopping_pref::Vector{Float64} = []
        ms_hop_term_i_jprime = MajoranaString(2 * N_spinful_sites, [4 * i - 3, 4 * j - 2])
        ms_hop_term_iprime_j = MajoranaString(2 * N_spinful_sites, [4 * i - 2, 4 * j - 3])

        push!(hopping_ms, ms_hop_term_i_jprime)
        push!(hopping_ms, ms_hop_term_iprime_j)
        push!(hopping_pref, +1)
        push!(hopping_pref, -1)

        push!(circ_up_hopping, FermionicGate(hopping_ms, hopping_pref))
        push!(thetas_up_hopping, -t * dt)

        if return_mps_instructions
            push!(mps_instructions, ("h", 2 * i - 1, 2 * j - 1))
            push!(mps_thetas, -dt * t)
        end
    end

    #repulsion term
    circ_repulsion::Vector{FermionicGate} = []
    thetas_repulsion::Vector{Float64} = []
    for k = 1:N_spinful_sites
        repulsion_ms::Vector{MajoranaString} = []
        repulsion_pref::Vector{Float64} = []
        ms_num_term_up = MajoranaString(2 * N_spinful_sites, [4 * k - 3, 4 * k - 2])
        ms_num_term_down = MajoranaString(2 * N_spinful_sites, [4 * k - 1, 4 * k])
        #@show ms_num_term_up, ms_num_term_down
        up_down_pref, ms_num_term_up_down = ms_mult(ms_num_term_up, ms_num_term_down)
        #@show up_down_pref, ms_num_term_up_down
        @assert imag(up_down_pref) < 1.e-12
        up_down_pref = real(up_down_pref)

        push!(repulsion_ms, ms_num_term_up)
        push!(repulsion_ms, ms_num_term_down)
        push!(repulsion_ms, ms_num_term_up_down)
        push!(repulsion_pref, 1.)
        push!(repulsion_pref, 1.)
        push!(repulsion_pref, up_down_pref)

        push!(circ_repulsion, FermionicGate(repulsion_ms, repulsion_pref))
        push!(thetas_repulsion, U * dt / 2)

        if return_mps_instructions
            push!(mps_instructions, ("nund", 2 * k - 1, 2 * k))
            push!(mps_thetas, dt * U)
        end

    end

    if return_separated
        circs = [circ_down_hopping, circ_up_hopping, circ_repulsion]
        thetas = [thetas_down_hopping, thetas_up_hopping, thetas_repulsion]
    else
        circs = vcat(circ_down_hopping, circ_up_hopping, circ_repulsion)
        thetas = vcat(thetas_down_hopping, thetas_up_hopping, thetas_repulsion)
    end

    if return_mps_instructions
        return circs, thetas, mps_instructions, mps_thetas
    end
    return circs, thetas
end


function hubbard_circ_fermionic_sites(topology, N_spinful_sites::Int, n_layers::Int, t::Float64, U::Float64, T::Float64; return_mps_instructions=false)
    circ::Vector{MajoranaRotation} = []
    thetas::Vector{Float64} = []
    mps_instructions = []
    mps_thetas = []
    dt = T / n_layers

    if return_mps_instructions
        circ_single_layer, thetas_single_layer, mps_instructions_single_layer, mps_thetas_single_layer = hubbard_circ_fermionic_sites_single_layer(topology, N_spinful_sites, t, U, dt; return_mps_instructions=true, return_separated=false)
        circ = repeat(circ_single_layer, n_layers)
        thetas = repeat(thetas_single_layer, n_layers)
        mps_instructions = repeat(mps_instructions_single_layer, n_layers)
        mps_thetas = repeat(mps_thetas_single_layer, n_layers)
        return circ, thetas, mps_instructions, mps_thetas
    end

    circ_single_layer, thetas_single_layer = hubbard_circ_fermionic_sites_single_layer(topology, N_spinful_sites, t, U, dt; return_mps_instructions=false, return_separated=false)
    circ = repeat(circ_single_layer, n_layers)
    thetas = repeat(thetas_single_layer, n_layers)

    return circ, thetas
end

function hubbard_circ_fermionic_sites_second_order_single_layer(topology, N_spinful_sites::Int, t::Float64, U::Float64, dt::Float64; return_mps_instructions=false, return_separated=false)
    mps_instructions = []
    mps_thetas = []

    #down part 
    circ_down_hopping_initial::Vector{MajoranaRotation} = []
    thetas_down_hopping_initial::Vector{Float64} = []
    for (i, j) in topology
        @assert i < j
        #@show i, j
        ms_hop_term_i_jprime = MajoranaString(2 * N_spinful_sites, [4 * i - 1, 4 * j])
        ms_hop_term_iprime_j = MajoranaString(2 * N_spinful_sites, [4 * i, 4 * j - 1])

        push!(circ_down_hopping_initial, MajoranaRotation(ms_hop_term_i_jprime))
        push!(thetas_down_hopping_initial, -t * dt / 2)
        push!(circ_down_hopping_initial, MajoranaRotation(ms_hop_term_iprime_j))
        push!(thetas_down_hopping_initial, +t * dt / 2)

        if return_mps_instructions
            #push!(mps_instructions, ("h", 2 * i, 2 * j))
            push!(mps_instructions, ("hd", i, j))
            push!(mps_thetas, -dt * t / 2)
        end
    end

    #up part 
    circ_up_hopping_initial::Vector{MajoranaRotation} = []
    thetas_up_hopping_initial::Vector{Float64} = []
    for (i, j) in topology
        @assert i < j
        ms_hop_term_i_jprime = MajoranaString(2 * N_spinful_sites, [4 * i - 3, 4 * j - 2])
        ms_hop_term_iprime_j = MajoranaString(2 * N_spinful_sites, [4 * i - 2, 4 * j - 3])

        push!(circ_up_hopping_initial, MajoranaRotation(ms_hop_term_i_jprime))
        push!(thetas_up_hopping_initial, -t * dt / 2)
        push!(circ_up_hopping_initial, MajoranaRotation(ms_hop_term_iprime_j))
        push!(thetas_up_hopping_initial, +t * dt / 2)
        if return_mps_instructions
            #push!(mps_instructions, ("h", 2 * i - 1, 2 * j - 1))
            push!(mps_instructions, ("hu", i, j))
            push!(mps_thetas, -dt * t / 2)
        end
    end

    #repulsion term 
    circ_repulsion::Vector{MajoranaRotation} = []
    thetas_repulsion::Vector{Float64} = []
    for k = 1:N_spinful_sites
        ms_num_term_up = MajoranaString(2 * N_spinful_sites, [4 * k - 3, 4 * k - 2])
        ms_num_term_down = MajoranaString(2 * N_spinful_sites, [4 * k - 1, 4 * k])
        #@show ms_num_term_up, ms_num_term_down
        up_down_pref, ms_num_term_up_down = ms_mult(ms_num_term_up, ms_num_term_down)
        #@show up_down_pref, ms_num_term_up_down
        @assert imag(up_down_pref) < 1.e-12
        up_down_pref = real(up_down_pref)

        push!(circ_repulsion, MajoranaRotation(ms_num_term_up))
        push!(thetas_repulsion, U * dt / 2)
        push!(circ_repulsion, MajoranaRotation(ms_num_term_down))
        push!(thetas_repulsion, U * dt / 2)
        push!(circ_repulsion, MajoranaRotation(ms_num_term_up_down))
        push!(thetas_repulsion, up_down_pref * U * dt / 2)

        if return_mps_instructions
            #push!(mps_instructions, ("nund", 2 * k - 1, 2 * k))
            push!(mps_instructions, ("nund", k))
            push!(mps_thetas, dt * U)
        end

    end

    #up part
    circ_up_hopping_final::Vector{MajoranaRotation} = []
    thetas_up_hopping_final::Vector{Float64} = []
    for (i, j) in reverse(topology)
        @assert i < j
        ms_hop_term_i_jprime = MajoranaString(2 * N_spinful_sites, [4 * i - 3, 4 * j - 2])
        ms_hop_term_iprime_j = MajoranaString(2 * N_spinful_sites, [4 * i - 2, 4 * j - 3])

        push!(circ_up_hopping_final, MajoranaRotation(ms_hop_term_iprime_j))
        push!(thetas_up_hopping_final, +t * dt / 2)
        push!(circ_up_hopping_final, MajoranaRotation(ms_hop_term_i_jprime))
        push!(thetas_up_hopping_final, -t * dt / 2)

        if return_mps_instructions
            #push!(mps_instructions, ("h", 2 * i - 1, 2 * j - 1))
            push!(mps_instructions, ("hu", i, j))
            push!(mps_thetas, -dt * t / 2)
        end
    end

    #down part
    circ_down_hopping_final::Vector{MajoranaRotation} = []
    thetas_down_hopping_final::Vector{Float64} = []
    for (i, j) in reverse(topology)
        @assert i < j
        #@show i, j
        ms_hop_term_i_jprime = MajoranaString(2 * N_spinful_sites, [4 * i - 1, 4 * j])
        ms_hop_term_iprime_j = MajoranaString(2 * N_spinful_sites, [4 * i, 4 * j - 1])

        push!(circ_down_hopping_final, MajoranaRotation(ms_hop_term_iprime_j))
        push!(thetas_down_hopping_final, +t * dt / 2)
        push!(circ_down_hopping_final, MajoranaRotation(ms_hop_term_i_jprime))
        push!(thetas_down_hopping_final, -t * dt / 2)
        if return_mps_instructions
            #push!(mps_instructions, ("h", 2 * i, 2 * j))
            push!(mps_instructions, ("hd", i, j))
            push!(mps_thetas, -dt * t / 2)
        end
    end

    if return_separated
        circs = [circ_down_hopping_initial, circ_up_hopping_initial, circ_repulsion, circ_up_hopping_final, circ_down_hopping_final]
        thetas = [thetas_down_hopping_initial, thetas_up_hopping_initial, thetas_repulsion, thetas_up_hopping_final, thetas_down_hopping_final]
    else
        circs = vcat(circ_down_hopping_initial, circ_up_hopping_initial, circ_repulsion, circ_up_hopping_final, circ_down_hopping_final)
        thetas = vcat(thetas_down_hopping_initial, thetas_up_hopping_initial, thetas_repulsion, thetas_up_hopping_final, thetas_down_hopping_final)
    end

    if return_mps_instructions
        return circs, thetas, mps_instructions, mps_thetas
    end
    return circs, thetas

end

function fermionic_hubbard_circ_fermionic_sites_second_order_single_layer(topology, N_spinful_sites::Int, t::Float64, U::Float64, dt::Float64; return_mps_instructions=false, return_separated=false)
    mps_instructions = []
    mps_thetas = []

    #down part 
    circ_down_hopping_initial::Vector{FermionicGate} = []
    thetas_down_hopping_initial::Vector{Float64} = []
    for (i, j) in topology
        @assert i < j
        #@show i, j
        hopping_ms::Vector{MajoranaString} = []
        hopping_pref::Vector{Float64} = []
        ms_hop_term_i_jprime = MajoranaString(2 * N_spinful_sites, [4 * i - 1, 4 * j])
        ms_hop_term_iprime_j = MajoranaString(2 * N_spinful_sites, [4 * i, 4 * j - 1])

        push!(hopping_ms, ms_hop_term_i_jprime)
        push!(hopping_ms, ms_hop_term_iprime_j)
        push!(hopping_pref, +1)
        push!(hopping_pref, -1)

        push!(circ_down_hopping_initial, FermionicGate(hopping_ms, hopping_pref))
        push!(thetas_down_hopping_initial, -t * dt)

        if return_mps_instructions
            push!(mps_instructions, ("h", 2 * i, 2 * j))
            push!(mps_thetas, -dt * t / 2)
        end
    end

    #up part 
    circ_up_hopping_initial::Vector{FermionicGate} = []
    thetas_up_hopping_initial::Vector{Float64} = []
    for (i, j) in topology
        @assert i < j
        hopping_ms::Vector{MajoranaString} = []
        hopping_pref::Vector{Float64} = []
        ms_hop_term_i_jprime = MajoranaString(2 * N_spinful_sites, [4 * i - 3, 4 * j - 2])
        ms_hop_term_iprime_j = MajoranaString(2 * N_spinful_sites, [4 * i - 2, 4 * j - 3])

        push!(hopping_ms, ms_hop_term_i_jprime)
        push!(hopping_ms, ms_hop_term_iprime_j)
        push!(hopping_pref, +1)
        push!(hopping_pref, -1)

        push!(circ_up_hopping_initial, FermionicGate(hopping_ms, hopping_pref))
        push!(thetas_up_hopping_initial, -t * dt)

        if return_mps_instructions
            push!(mps_instructions, ("h", 2 * i - 1, 2 * j - 1))
            push!(mps_thetas, -dt * t / 2)
        end
    end

    #repulsion term 
    circ_repulsion::Vector{FermionicGate} = []
    thetas_repulsion::Vector{Float64} = []
    for k = 1:N_spinful_sites
        repulsion_ms::Vector{MajoranaString} = []
        repulsion_pref::Vector{Float64} = []
        ms_num_term_up = MajoranaString(2 * N_spinful_sites, [4 * k - 3, 4 * k - 2])
        ms_num_term_down = MajoranaString(2 * N_spinful_sites, [4 * k - 1, 4 * k])
        #@show ms_num_term_up, ms_num_term_down
        up_down_pref, ms_num_term_up_down = ms_mult(ms_num_term_up, ms_num_term_down)
        #@show up_down_pref, ms_num_term_up_down
        @assert imag(up_down_pref) < 1.e-12
        up_down_pref = real(up_down_pref)

        push!(repulsion_ms, ms_num_term_up)
        push!(repulsion_ms, ms_num_term_down)
        push!(repulsion_ms, ms_num_term_up_down)
        push!(repulsion_pref, 1.)
        push!(repulsion_pref, 1.)
        push!(repulsion_pref, up_down_pref)

        push!(circ_repulsion, FermionicGate(repulsion_ms, repulsion_pref))
        push!(thetas_repulsion, U * dt / 2)

        if return_mps_instructions
            push!(mps_instructions, ("nund", 2 * k - 1, 2 * k))
            push!(mps_thetas, dt * U)
        end

    end

    #up part
    circ_up_hopping_final::Vector{FermionicGate} = []
    thetas_up_hopping_final::Vector{Float64} = []
    for (i, j) in reverse(topology)
        @assert i < j
        hopping_ms::Vector{MajoranaString} = []
        hopping_pref::Vector{Float64} = []
        ms_hop_term_i_jprime = MajoranaString(2 * N_spinful_sites, [4 * i - 3, 4 * j - 2])
        ms_hop_term_iprime_j = MajoranaString(2 * N_spinful_sites, [4 * i - 2, 4 * j - 3])

        push!(hopping_ms, ms_hop_term_i_jprime)
        push!(hopping_ms, ms_hop_term_iprime_j)
        push!(hopping_pref, +1)
        push!(hopping_pref, -1)

        push!(circ_up_hopping_final, FermionicGate(hopping_ms, hopping_pref))
        push!(thetas_up_hopping_final, -t * dt)

        if return_mps_instructions
            push!(mps_instructions, ("h", 2 * i - 1, 2 * j - 1))
            push!(mps_thetas, -dt * t / 2)
        end
    end

    #down part
    circ_down_hopping_final::Vector{FermionicGate} = []
    thetas_down_hopping_final::Vector{Float64} = []
    for (i, j) in reverse(topology)
        @assert i < j
        hopping_ms::Vector{MajoranaString} = []
        hopping_pref::Vector{Float64} = []
        ms_hop_term_i_jprime = MajoranaString(2 * N_spinful_sites, [4 * i - 1, 4 * j])
        ms_hop_term_iprime_j = MajoranaString(2 * N_spinful_sites, [4 * i, 4 * j - 1])

        push!(hopping_ms, ms_hop_term_i_jprime)
        push!(hopping_ms, ms_hop_term_iprime_j)
        push!(hopping_pref, +1)
        push!(hopping_pref, -1)

        push!(circ_down_hopping_final, FermionicGate(hopping_ms, hopping_pref))
        push!(thetas_down_hopping_final, -t * dt)

        if return_mps_instructions
            push!(mps_instructions, ("h", 2 * i, 2 * j))
            push!(mps_thetas, -dt * t / 2)
        end
    end

    if return_separated
        circs = [circ_down_hopping_initial, circ_up_hopping_initial, circ_repulsion, circ_up_hopping_final, circ_down_hopping_final]
        thetas = [thetas_down_hopping_initial, thetas_up_hopping_initial, thetas_repulsion, thetas_up_hopping_final, thetas_down_hopping_final]
    else
        circs = vcat(circ_down_hopping_initial, circ_up_hopping_initial, circ_repulsion, circ_up_hopping_final, circ_down_hopping_final)
        thetas = vcat(thetas_down_hopping_initial, thetas_up_hopping_initial, thetas_repulsion, thetas_up_hopping_final, thetas_down_hopping_final)
    end

    if return_mps_instructions
        return circs, thetas, mps_instructions, mps_thetas
    end
    return circs, thetas

end



function hubbard_circ_fermionic_sites_second_order(topology, N_spinful_sites::Int, n_layers::Int, t::Float64, U::Float64, T::Float64; return_mps_instructions=false)
    circ::Vector{MajoranaRotation} = []
    thetas::Vector{Float64} = []
    mps_instructions = []
    mps_thetas = []
    dt = T / n_layers

    if return_mps_instructions
        circ_single_layer, thetas_single_layer, mps_instructions_single_layer, mps_thetas_single_layer = hubbard_circ_fermionic_sites_second_order_single_layer(topology, N_spinful_sites, t, U, dt; return_mps_instructions=true, return_separated=false)
        circ = repeat(circ_single_layer, n_layers)
        thetas = repeat(thetas_single_layer, n_layers)
        mps_instructions = repeat(mps_instructions_single_layer, n_layers)
        mps_thetas = repeat(mps_thetas_single_layer, n_layers)
        return circ, thetas, mps_instructions, mps_thetas
    end

    circ_single_layer, thetas_single_layer = hubbard_circ_fermionic_sites_second_order_single_layer(topology, N_spinful_sites, t, U, dt; return_mps_instructions=false, return_separated=false)
    circ = repeat(circ_single_layer, n_layers)
    thetas = repeat(thetas_single_layer, n_layers)

    return circ, thetas
end