"""
    hubbard_circ_fermionic_sites_single_layer(topology, N_spinful_sites::Int, t::Float64, U::Float64, dt::Float64; return_mps_instructions=false, return_separated=false)

Build one first-order Trotter layer of the time evolution ``e^{-i \\, dt \\, H}`` of the Fermi-Hubbard Hamiltonian ``H = -t \\sum_{(i,j), \\sigma} (f_{i\\sigma}^\\dagger f_{j\\sigma} + f_{j\\sigma}^\\dagger f_{i\\sigma}) + U \\sum_i n_{i\\uparrow} n_{i\\downarrow}`` as a circuit of `MajoranaRotation`s.
The bonds ``(i, j)`` with ``i < j`` are given by `topology`, on a system of `N_spinful_sites` spinful sites.
Returns a tuple `(circuit, thetas)`, with the layer ordered as spin-down hoppings, spin-up hoppings, on-site repulsion.
If `return_separated=true`, `circuit` and `thetas` are instead grouped into these three blocks as vectors of vectors.
If `return_mps_instructions=true`, additionally returns `(mps_instructions, mps_thetas)` describing the equivalent gate sequence for an MPS simulation.
"""
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

"""
    fermionic_hubbard_circ_fermionic_sites_single_layer(topology, N_spinful_sites::Int, t::Float64, U::Float64, dt::Float64; return_mps_instructions=false, return_separated=false)

Variant of `hubbard_circ_fermionic_sites_single_layer` that groups the Majorana strings of each Hamiltonian term into `FermionicRotation` gates instead of individual `MajoranaRotation`s.
Same arguments and return values.
"""
function fermionic_hubbard_circ_fermionic_sites_single_layer(topology, N_spinful_sites::Int, t::Float64, U::Float64, dt::Float64; return_mps_instructions=false, return_separated=false)
    mps_instructions = []
    mps_thetas = []

    #down part
    circ_down_hopping::Vector{FermionicRotation} = []
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

        push!(circ_down_hopping, FermionicRotation(hopping_ms, hopping_pref))
        push!(thetas_down_hopping, -t * dt)

        if return_mps_instructions
            push!(mps_instructions, ("h", 2 * i, 2 * j))
            push!(mps_thetas, -dt * t)
        end
    end

    #up part
    circ_up_hopping::Vector{FermionicRotation} = []
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

        push!(circ_up_hopping, FermionicRotation(hopping_ms, hopping_pref))
        push!(thetas_up_hopping, -t * dt)

        if return_mps_instructions
            push!(mps_instructions, ("h", 2 * i - 1, 2 * j - 1))
            push!(mps_thetas, -dt * t)
        end
    end

    #repulsion term
    circ_repulsion::Vector{FermionicRotation} = []
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

        push!(circ_repulsion, FermionicRotation(repulsion_ms, repulsion_pref))
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


"""
    hubbard_circ_fermionic_sites(topology, N_spinful_sites::Int, n_layers::Int, t::Float64, U::Float64, T::Float64; return_mps_instructions=false)

Build the first-order Trotter circuit for the Fermi-Hubbard time evolution ``e^{-i \\, T \\, H}`` by repeating `hubbard_circ_fermionic_sites_single_layer` with `dt = T / n_layers` for `n_layers` layers.
Returns `(circuit, thetas)`, and additionally `(mps_instructions, mps_thetas)` if `return_mps_instructions=true`.
"""
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

"""
    hubbard_circ_fermionic_sites_second_order_single_layer(topology, N_spinful_sites::Int, t::Float64, U::Float64, dt::Float64; return_mps_instructions=false, return_separated=false)

Build one second-order (symmetric) Trotter layer of the Fermi-Hubbard time evolution ``e^{-i \\, dt \\, H}`` as a circuit of `MajoranaRotation`s: half-step spin-down and spin-up hoppings, a full repulsion step, then the up and down hoppings again in reversed order for the second half step.
Same arguments and return values as `hubbard_circ_fermionic_sites_single_layer`; with `return_separated=true` the circuit is grouped into these five blocks.
"""
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

"""
    fermionic_hubbard_circ_fermionic_sites_second_order_single_layer(topology, N_spinful_sites::Int, t::Float64, U::Float64, dt::Float64; return_mps_instructions=false, return_separated=false)

Variant of `hubbard_circ_fermionic_sites_second_order_single_layer` that groups the Majorana strings of each Hamiltonian term into `FermionicRotation` gates instead of individual `MajoranaRotation`s.
Same arguments and return values.
"""
function fermionic_hubbard_circ_fermionic_sites_second_order_single_layer(topology, N_spinful_sites::Int, t::Float64, U::Float64, dt::Float64; return_mps_instructions=false, return_separated=false)
    mps_instructions = []
    mps_thetas = []

    #down part 
    circ_down_hopping_initial::Vector{FermionicRotation} = []
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

        push!(circ_down_hopping_initial, FermionicRotation(hopping_ms, hopping_pref))
        push!(thetas_down_hopping_initial, -t * dt)

        if return_mps_instructions
            push!(mps_instructions, ("h", 2 * i, 2 * j))
            push!(mps_thetas, -dt * t / 2)
        end
    end

    #up part 
    circ_up_hopping_initial::Vector{FermionicRotation} = []
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

        push!(circ_up_hopping_initial, FermionicRotation(hopping_ms, hopping_pref))
        push!(thetas_up_hopping_initial, -t * dt)

        if return_mps_instructions
            push!(mps_instructions, ("h", 2 * i - 1, 2 * j - 1))
            push!(mps_thetas, -dt * t / 2)
        end
    end

    #repulsion term 
    circ_repulsion::Vector{FermionicRotation} = []
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

        push!(circ_repulsion, FermionicRotation(repulsion_ms, repulsion_pref))
        push!(thetas_repulsion, U * dt / 2)

        if return_mps_instructions
            push!(mps_instructions, ("nund", 2 * k - 1, 2 * k))
            push!(mps_thetas, dt * U)
        end

    end

    #up part
    circ_up_hopping_final::Vector{FermionicRotation} = []
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

        push!(circ_up_hopping_final, FermionicRotation(hopping_ms, hopping_pref))
        push!(thetas_up_hopping_final, -t * dt)

        if return_mps_instructions
            push!(mps_instructions, ("h", 2 * i - 1, 2 * j - 1))
            push!(mps_thetas, -dt * t / 2)
        end
    end

    #down part
    circ_down_hopping_final::Vector{FermionicRotation} = []
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

        push!(circ_down_hopping_final, FermionicRotation(hopping_ms, hopping_pref))
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



"""
    hubbard_circ_fermionic_sites_second_order(topology, N_spinful_sites::Int, n_layers::Int, t::Float64, U::Float64, T::Float64; return_mps_instructions=false)

Build the second-order Trotter circuit for the Fermi-Hubbard time evolution ``e^{-i \\, T \\, H}`` by repeating `hubbard_circ_fermionic_sites_second_order_single_layer` with `dt = T / n_layers` for `n_layers` layers.
Returns `(circuit, thetas)`, and additionally `(mps_instructions, mps_thetas)` if `return_mps_instructions=true`.
"""
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