using MajoranaPropagation
using PauliPropagation
using ITensors
using ITensorMPS

using Random
Random.seed!(42)

function string_to_mpo_spinless(ms::TT, n_fermions, sites) where {TT<:Integer}
    ms_mpo = MPO(sites)

    all_operators = Int[]

    for i = 1:n_fermions
        gamma = ((ms >> (2 * i - 2)) & TT(1))
        gamma_prime = ((ms >> (2 * i - 1)) & TT(1))
        if gamma == 1
            push!(all_operators, i)
        end
        if gamma_prime == 1
            push!(all_operators, -i)
        end
    end

    for k = 1:2:length(all_operators)
        o1 = all_operators[k]
        o2 = all_operators[k+1]
        os = OpSum()
        if o1 + o2 == 0
            i = abs(o1)
            os += -2im, "N", i
            os += 1im, "I", i
        elseif o1 > 0
            if o2 > 0
                # gamma_i gamma_j
                i = abs(o1)
                j = abs(o2)
                os += "Cdag", i, "Cdag", j
                os += "Cdag", i, "C", j
                os += "C", i, "Cdag", j
                os += "C", i, "C", j
            else
                # gamma_i gamma'_j
                i = abs(o1)
                j = abs(o2)
                os += 1im, "Cdag", i, "Cdag", j
                os += -1im, "Cdag", i, "C", j
                os += 1im, "C", i, "Cdag", j
                os += -1im, "C", i, "C", j
            end
        else
            if o2 > 0
                # gamma'_i gamma_j
                i = abs(o1)
                j = abs(o2)
                os += 1im, "Cdag", i, "Cdag", j
                os += 1im, "Cdag", i, "C", j
                os += -1im, "C", i, "Cdag", j
                os += -1im, "C", i, "C", j
            else
                # gamma'_i gamma'_j
                i = abs(o1)
                j = abs(o2)
                os += -1, "Cdag", i, "Cdag", j
                os += 1, "Cdag", i, "C", j
                os += 1, "C", i, "Cdag", j
                os += -1, "C", i, "C", j
            end
        end
        os_mpo = MPO(os, sites)
        #@show os 
        if k == 1
            ms_mpo = os_mpo
        else
            ms_mpo = ITensorMPS.apply(ms_mpo, os_mpo)
        end
    end
    ms_mpo *= (1im)^(omega_L_mult(ms))
    return ms_mpo
end


function string_to_mpo_spinful(ms::TT, n_sites, sites) where {TT<:Integer}
    ms_mpo = MPO(sites)

    all_operators = []

    for i = 1:n_sites
        gamma_up = ((ms >> (4 * i - 4)) & TT(1))
        gamma_prime_up = ((ms >> (4 * i - 3)) & TT(1))
        gamma_dn = ((ms >> (4 * i - 2)) & TT(1))
        gamma_prime_dn = ((ms >> (4 * i - 1)) & TT(1))
        if gamma_up == 1
            push!(all_operators, (i, "up"))
        end
        if gamma_prime_up == 1
            push!(all_operators, (-i, "up"))
        end
        if gamma_dn == 1
            push!(all_operators, (i, "dn"))
        end
        if gamma_prime_dn == 1
            push!(all_operators, (-i, "dn"))
        end
    end

    for k = 1:2:length(all_operators)
        o1, sigma1 = all_operators[k]
        o2, sigma2 = all_operators[k+1]
        os = OpSum()
        if o1 + o2 == 0 && sigma1 == sigma2
            i = abs(o1)
            os += -2im, "N$sigma1", i
            os += 1im, "I", i
        elseif o1 > 0
            if o2 > 0
                # gamma_i gamma_j
                i = abs(o1)
                j = abs(o2)
                os += "Cdag$sigma1", i, "Cdag$sigma2", j
                os += "Cdag$sigma1", i, "C$sigma2", j
                os += "C$sigma1", i, "Cdag$sigma2", j
                os += "C$sigma1", i, "C$sigma2", j
            else
                # gamma_i gamma'_j
                i = abs(o1)
                j = abs(o2)
                os += 1im, "Cdag$sigma1", i, "Cdag$sigma2", j
                os += -1im, "Cdag$sigma1", i, "C$sigma2", j
                os += 1im, "C$sigma1", i, "Cdag$sigma2", j
                os += -1im, "C$sigma1", i, "C$sigma2", j
            end
        else
            if o2 > 0
                # gamma'_i gamma_j
                i = abs(o1)
                j = abs(o2)
                os += 1im, "Cdag$sigma1", i, "Cdag$sigma2", j
                os += 1im, "Cdag$sigma1", i, "C$sigma2", j
                os += -1im, "C$sigma1", i, "Cdag$sigma2", j
                os += -1im, "C$sigma1", i, "C$sigma2", j
            else
                # gamma'_i gamma'_j
                i = abs(o1)
                j = abs(o2)
                os += -1, "Cdag$sigma1", i, "Cdag$sigma2", j
                os += 1, "Cdag$sigma1", i, "C$sigma2", j
                os += 1, "C$sigma1", i, "Cdag$sigma2", j
                os += -1, "C$sigma1", i, "C$sigma2", j
            end
        end
        os_mpo = MPO(os, sites)
        #@show os 
        if k == 1
            ms_mpo = os_mpo
        else
            ms_mpo = ITensorMPS.apply(ms_mpo, os_mpo)
        end
    end
    ms_mpo *= (1im)^(omega_L_mult(ms))
    return ms_mpo
end

function make_mps_sites(n_sites, occupied_sites_up, occupied_sites_dn)
    mps_sites = []
    for i in 1:n_sites
        if i in occupied_sites_up && i in occupied_sites_dn
            push!(mps_sites, "UpDn")
        elseif i in occupied_sites_up
            push!(mps_sites, "Up")
        elseif i in occupied_sites_dn
            push!(mps_sites, "Dn")
        else
            push!(mps_sites, "Emp")
        end
    end
    return mps_sites
end

function spinless()
    n_fermions = 6
    sites = siteinds("Fermion", n_fermions)

    max_elements_in_superposition = 20
    max_particles = min(10, n_fermions)
    n_tests = 100
    n_strings_per_test = 200

    TT = getinttype(n_fermions)
    max_val = MajoranaString(n_fermions, [2 * n_fermions - 1]).gammas

    non_zero_overlaps = 0

    for l in 1:n_tests
        n_focks = rand(1:max_elements_in_superposition)
        n_particles = rand(1:max_particles)
        fock_states = Vector{FockState}(undef, n_focks)
        fock_states_mps = Vector{MPS}(undef, n_focks)

        occupied_sites_orig = randperm(n_fermions)[1:n_particles-1]

        for j in 1:n_focks
            occupied_sites = deepcopy(occupied_sites_orig)
            while length(occupied_sites) < n_particles
                rand_site = rand(1:n_fermions)
                if !(rand_site in occupied_sites)
                    push!(occupied_sites, rand_site)
                end
            end
            occupied_sites = sort(occupied_sites)
            fock_states[j] = FockState(n_fermions, occupied_sites)
            fock_states_mps[j] = MPS(sites, [i in occupied_sites ? "1" : "0" for i in 1:n_fermions])
            #@show fock_states[j]
        end


        # check overlaps of the form <fock1 | ms | fock2> for random Majorana strings and random Fock states
        k = 1
        while k <= n_strings_per_test
            rand_ms = rand(TT)
            if rand_ms > max_val || get_weight(rand_ms) % 2 != 0 || rand_ms == 0
                continue
            end
            f1_index = rand(1:n_focks)
            f2_index = rand(1:n_focks)
            f1 = fock_states[f1_index]
            f2 = fock_states[f2_index]
            ms_mpo = string_to_mpo_spinless(rand_ms, n_fermions, sites)
            MP_overlap = overlapwithfock(rand_ms, f1, f2, n_fermions)
            mps_overlap = inner(fock_states_mps[f1_index]', ms_mpo, fock_states_mps[f2_index])
            #@show MP_overlap, mps_overlap, bitstring(rand_ms)
            @assert abs(MP_overlap - mps_overlap) < 1.e-12
            if abs(MP_overlap) > 1.e-12
                #@show MP_overlap, mps_overlap
                non_zero_overlaps += 1
            end
            k += 1
        end

        #println("-----")

        #check overlaps of superpositions 
        k = 1
        while k <= n_strings_per_test
            rand_ms = rand(TT)
            if rand_ms > max_val || get_weight(rand_ms) % 2 != 0 || rand_ms == 0
                continue
            end
            msum = MajoranaSum(Float64, n_fermions, false)
            set!(msum, rand_ms, 1.)
            superposition_coefficients = randn(ComplexF64, n_focks)
            superposition_coefficients ./= sqrt(sum(abs2, superposition_coefficients)) # normalize
            ms_mpo = string_to_mpo_spinless(rand_ms, n_fermions, sites)
            MP_overlap = overlapwithfock(msum, fock_states, superposition_coefficients)
            psi_superposition = sum(superposition_coefficients[j] * fock_states_mps[j] for j in 1:n_focks)
            mps_overlap = inner(psi_superposition', ms_mpo, psi_superposition)
            @assert abs(MP_overlap - mps_overlap) < 1.e-12
            if abs(MP_overlap) > 1.e-12
                #@show MP_overlap, mps_overlap
                non_zero_overlaps += 1
            end
            k += 1
        end
    end
    println("All tests passed! Non-zero overlaps computed: ", non_zero_overlaps)
end

function spinful()
    n_sites = 5
    n_fermions = 2 * n_sites
    sites = siteinds("Electron", n_sites)

    max_elements_in_superposition = 20
    max_particles = min(10, n_sites)
    n_tests = 100
    n_strings_per_test = 200

    TT = getinttype(n_fermions)
    max_val = MajoranaString(n_fermions, [2 * n_fermions - 1]).gammas

    non_zero_overlaps = 0

    for l in 1:n_tests
        n_focks = rand(1:max_elements_in_superposition)
        n_particles_up = rand(1:max_particles)
        n_particles_dn = rand(1:max_particles)

        fock_states = Vector{FockState}(undef, n_focks)
        fock_states_mps = Vector{MPS}(undef, n_focks)

        occupied_sites_orig_up = randperm(n_sites)[1:n_particles_up-1]
        occupied_sites_orig_dn = randperm(n_sites)[1:n_particles_dn-1]

        for j in 1:n_focks
            occupied_sites_up = deepcopy(occupied_sites_orig_up)
            occupied_sites_dn = deepcopy(occupied_sites_orig_dn)
            while length(occupied_sites_up) < n_particles_up
                rand_site = rand(1:n_sites)
                if !(rand_site in occupied_sites_up)
                    push!(occupied_sites_up, rand_site)
                end
            end
            while length(occupied_sites_dn) < n_particles_dn
                rand_site = rand(1:n_sites)
                if !(rand_site in occupied_sites_dn)
                    push!(occupied_sites_dn, rand_site)
                end
            end

            fock_states[j] = FockState(n_sites, occupied_sites_up, occupied_sites_dn)
            mps_sites = make_mps_sites(n_sites, occupied_sites_up, occupied_sites_dn)

            fock_states_mps[j] = MPS(sites, mps_sites)
        end


        # check overlaps of the form <fock1 | ms | fock2> for random Majorana strings and random Fock states
        k = 1
        unpaired_mask = create_unpaired_mask(2 * n_sites)
        while k <= n_strings_per_test
            rand_ms = rand(TT)
            if rand_ms > max_val || get_weight(rand_ms) % 2 != 0 || rand_ms == 0
                continue
            end
            f1_index = rand(1:n_focks)
            f2_index = rand(1:n_focks)
            f1 = fock_states[f1_index]
            f2 = fock_states[f2_index]
            ms_mpo = string_to_mpo_spinful(rand_ms, n_sites, sites)
            
            MP_overlap = overlapwithfock(rand_ms, f1, f2, n_fermions)
            mps_overlap = inner(fock_states_mps[f1_index]', ms_mpo, fock_states_mps[f2_index])
            #MP_overlap = overlapwithfock(rand_ms, f1, f1, n_fermions)
            #mps_overlap = inner(fock_states_mps[f1_index]', ms_mpo, fock_states_mps[f1_index])
            
            #@show MP_overlap, mps_overlap, bitstring(rand_ms), f1
            #@show overlapwithfock(rand_ms, unpaired_mask, f1)
            @assert abs(MP_overlap - mps_overlap) < 1.e-12
            if abs(MP_overlap) > 1.e-12
                #@show MP_overlap, mps_overlap
                non_zero_overlaps += 1
            end
            k += 1
        end

        continue

        #println("-----")

        #check overlaps of superpositions 
        k = 1
        while k <= n_strings_per_test
            rand_ms = rand(TT)
            if rand_ms > max_val || get_weight(rand_ms) % 2 != 0 || rand_ms == 0
                continue
            end
            msum = MajoranaSum(Float64, n_sites, true)
            set!(msum, rand_ms, 1.)
            superposition_coefficients = randn(ComplexF64, n_focks)
            superposition_coefficients ./= sqrt(sum(abs2, superposition_coefficients)) # normalize
            ms_mpo = string_to_mpo_spinful(rand_ms, n_sites, sites)
            MP_overlap = overlapwithfock(msum, fock_states, superposition_coefficients)
            psi_superposition = sum(superposition_coefficients[j] * fock_states_mps[j] for j in 1:n_focks)
            mps_overlap = inner(psi_superposition', ms_mpo, psi_superposition)
            @assert abs(MP_overlap - mps_overlap) < 1.e-12
            if abs(MP_overlap) > 1.e-12
                #@show MP_overlap, mps_overlap
                non_zero_overlaps += 1
            end
            k += 1
        end
    end
    println("All tests passed! Non-zero overlaps computed: ", non_zero_overlaps)

end

function simple_spinless()
    n_fermions = 10
    sites = siteinds("Fermion", n_fermions)
    msum = MajoranaSum(n_fermions, :nn, [1, 2])
    @show msum

    f1_l = collect(1:2:n_fermions)
    f2_l = deepcopy(f1_l)
    f2_l[1] = 2

    f1 = FockState(n_fermions, f1_l)
    f2 = FockState(n_fermions, f2_l)

    f1_mps = MPS(sites, [i in f1_l ? "1" : "0" for i in 1:n_fermions])
    f2_mps = MPS(sites, [i in f2_l ? "1" : "0" for i in 1:n_fermions])



    for (ms, coeff) in msum
        if ms > 0
            println("---")
            ms_mpo = string_to_mpo_spinless(ms, n_fermions, sites)

            MP_overlap = overlapwithfock(ms, f1, f2, n_fermions)
            mps_overlap = inner(f1_mps', ms_mpo, f2_mps)
            @show MP_overlap
            @show mps_overlap
        end
    end
end

function simple_spinful()
    n_sites = 10
    sites = siteinds("Electron", n_sites)
    msum = MajoranaSum(n_sites, :nupndn, 1)
    @show msum

    f1_l_up = collect(1:2:n_sites)
    f1_l_dn = collect(2:2:n_sites)

    f2_l = deepcopy(f1_l_up)
    f2_l[1] = 2

    f1 = FockState(n_sites, f1_l_up, f1_l_dn)
    f2 = FockState(n_sites, f2_l, [])

    f1_mps = MPS(sites, make_mps_sites(n_sites, f1_l_up, f1_l_dn))
    f2_mps = MPS(sites, make_mps_sites(n_sites, f2_l, []))



    for (ms, coeff) in msum
        if ms > 0
            println("---")
            @show bitstring(ms)
            ms_mpo = string_to_mpo_spinful(ms, n_sites, sites)

            MP_overlap = overlapwithfock(ms, f1, f1, n_sites)
            mps_overlap = inner(f1_mps', ms_mpo, f1_mps)
            @show MP_overlap
            @show mps_overlap
        end
    end
end

#simple_spinless()
#simple_spinful()

spinless()
spinful()