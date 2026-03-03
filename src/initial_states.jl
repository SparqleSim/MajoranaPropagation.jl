"""
    FockState
A struct to represent a Fock basis state.
"""
struct FockState{TT<:Integer}
    n_sites::Int
    is_spinful::Bool
    occupied_sites::TT
end

"""
    FockState(n_sites::Int, occupied_sites_iter)
Create a spinless Fock basis state given an iterable with the occupied sites.
"""
function FockState(n_sites::Int, occupied_sites_iter::AbstractVector)
    TT = getinttype(n_sites)
    occupied_sites = _bitonesat(TT, (2 * site - 1 for site in occupied_sites_iter))
    return FockState(n_sites, false, occupied_sites)
end

"""
    FockState(n_sites::Int, up_occupied_sites_iter, down_occupied_sites_iter)
Create a spinful Fock basis state given an iterable with the occupied sites for spin-up and a list of occupied sites for spin-down fermions.
"""
function FockState(n_sites::Int, up_occupied_sites_iter::AbstractVector, down_occupied_sites_iter::AbstractVector)
    TT = getinttype(2 * n_sites)
    occupied_sites_list::Vector{Int} = []
    for site in up_occupied_sites_iter
        push!(occupied_sites_list, 2 * site - 1)
    end
    for site in down_occupied_sites_iter
        push!(occupied_sites_list, 2 * site)
    end
    sort!(occupied_sites_list)
    occupied_sites = _bitonesat(TT, (2 * site - 1 for site in occupied_sites_list))
    return FockState(n_sites, true, occupied_sites)
end

function FockState(n_sites::Integer, symb::Symbol, is_spinful::Bool; hole_positions=nothing, kwargs...)
    return FockState(n_sites, Val(symb), is_spinful; hole_positions=hole_positions, kwargs...)
end

function FockState(n_sites::Integer, ::Val{:checkerboard}, is_spinful::Bool; hole_positions=nothing, nx=-1, kwargs...)
    holes = isnothing(hole_positions) ? [] : hole_positions
    holes = isa(holes, Integer) ? [holes] : holes
    nx = nx == -1 ? n_sites : nx
    if is_spinful
        create_up_part_at::Vector{Int} = []
        create_down_part_at::Vector{Int} = []

        for j = 1:n_sites
            jy, jx = divrem(j - 1, nx)
            jy += 1
            jx += 1

            if j in holes
                continue
            elseif (jx % 2 == 1) && (jy % 2 == 1)
                push!(create_up_part_at, j)
            elseif (jx % 2 == 0) && (jy % 2 == 0)
                push!(create_up_part_at, j)
            else
                push!(create_down_part_at, j)
            end
        end
        return FockState(n_sites, create_up_part_at, create_down_part_at)
    else
        occupied_sites = [j for j in 1:n_sites if j % 2 == 1]
        occupied_sites = setdiff(occupied_sites, holes)
        return FockState(n_sites, occupied_sites)
    end
    @error "Invalid symbol for FockState constructor."
end

# printing for fock states
function Base.show(io::IO, fock_state::FockState)
    is_spinful = fock_state.is_spinful
    occupied_slots = fock_state.occupied_sites

    if is_spinful
        up_fermions::Vector{Integer} = []
        down_fermions::Vector{Integer} = []
        #loop over bits of occupied sites and separate up and down fermions
        for site in 1:fock_state.n_sites
            upbit = (occupied_slots >> (4 * site - 4)) & 1
            downbit = (occupied_slots >> (4 * site - 2)) & 1
            if upbit == 1
                push!(up_fermions, site)
            end
            if downbit == 1
                push!(down_fermions, site)
            end
        end

        print(io, "Fock state with $(length(up_fermions) + length(down_fermions)) fermions at positions\n    ↑: $(join(up_fermions, ", "))\n    ↓: $(join(down_fermions, ", "))\n")
    else
        occupied_sites::Vector{Integer} = []
        for site in 1:fock_state.n_sites
            bit = (occupied_slots >> (2 * site - 2)) & 1
            if bit == 1
                push!(occupied_sites, site)
            end
        end
        print(io, "Fock state with $(length(occupied_sites)) fermions at positions: $(join(occupied_sites, ", "))")
    end
end

"""
    fock_mask(msum::MajoranaSum)
Given a MajoranaSum, return a new MajoranaSum with terms that have an overlap with Fock states.
"""
function fock_mask(msum::MajoranaSum)
    clean_res = similar(msum)
    singles_filter = create_unpaired_mask(nfermions(msum))
    for (ms, coeff) in zip(terms(msum), coefficients(msum))
        if compute_unpaired(ms, singles_filter) > 0
            continue
        end
        set!(clean_res, ms, coeff)
    end
    return clean_res
end

"""
    overlapwithfock(msum::MajoranaSum, fock_state::FockState)
Compute the overlap <fock_state|msum|fock_state> where fock_state is a `FockState` object.
"""
function overlapwithfock(msum::AbstractMajoranaSum, fock_state::FockState)
    @assert is_spinful(msum) == fock_state.is_spinful "The MajoranaSum and the fock_state must both be spinful or both spinless."
    res = 0.
    unpaired_mask = create_unpaired_mask(nfermions(msum))
    for (ms, coeff) in zip(majoranas(msum), coefficients(msum))
        res += tonumber(coeff) * overlapwithfock(ms, unpaired_mask, fock_state)
    end
    return res
end

"""
    overlapwithfock(msum::MajoranaSum{TT,CT}, fock_state_1::FockState, fock_state_2::FockState) where {TT<:Integer,CT}

Evaluate the matrix element <fock_state_1|ms|fock_state_2> where fock_state_j are Fock basis states given as list of integers indicating which sites are occupied.
"""
function overlapwithfock(msum::MajoranaSum{TT,CT}, fock_state_1::FockState, fock_state_2::FockState) where {TT<:Integer,CT}
    @assert is_spinful(msum) == fock_state_1.is_spinful == fock_state_2.is_spinful "The MajoranaSum and the fock_states must both be spinful or both spinless."
    res = 0.
    n_fermions = nfermions(msum)
    for (ms, coeff) in zip(terms(msum), coefficients(msum))
        res += coeff * overlapwithfock(ms, fock_state_1, fock_state_2, n_fermions)
        #@show bitstring(ms), res
    end
    return res
end

"""
    overlapwithfock(ms::TT, unpaired_mask::TT, fock_state::FockState) where {TT<:Integer}
Compute the overlap <fock_state|ms|fock_state> where fock_state is a `FockState` object.
"""
function overlapwithfock(ms::TT, unpaired_mask::TT, fock_state::FockState) where {TT<:Integer}
    if compute_unpaired(ms, unpaired_mask) > 0
        return 0.
    end
    number_pref = get_weight(ms & fock_state.occupied_sites)
    sign = (1im)^omega_L_mult(ms) * (1im)^(get_weight(ms) / 2) * (-1)^number_pref
    return real(sign)
end

"""
    overlapwithfock(ms::TT, fock_state_1::FockState, fock_state_2::FockState) where {TT<:Integer}

Evaluate the matrix element <fock_state_1|ms|fock_state_2> where fock_state_j are Fock basis states given as list of integers indicating which sites are occupied.
"""
function overlapwithfock(ms::TT, fock_state_1::FockState, fock_state_2::FockState, n_fermions) where {TT<:Integer}
    res = (1im)^omega_L_mult(ms)
    for i = 1:n_fermions
        gamma = ((ms >> (2 * i - 2)) & TT(1))
        gamma_prime = ((ms >> (2 * i - 1)) & TT(1))
        if gamma == gamma_prime
            if (i in fock_state_2.occupied_sites) != (i in fock_state_1.occupied_sites)
                res *= 0.
                break
            else
                res *= (1im * (-1)^(i in fock_state_2.occupied_sites))^gamma
            end
        else
            if (i in fock_state_2.occupied_sites) == (i in fock_state_1.occupied_sites)
                res *= 0.
                break
            else
                res *= (1im * (-1)^(i in fock_state_2.occupied_sites))^gamma_prime * (-1)^(sum((j in fock_state_1.occupied_sites) for j = min(i + 1, n_fermions):n_fermions))
            end
        end
    end
    return res
end

"""
    overlapwithfock(msum::MajoranaSum, sites_with_particle_superposition::Vector{FockState}, superposition_coefficients::Vector{<:Union{Real,Complex}})
Compute the overlap <superposition|msum|superposition> where
- superposition is given as a vector of Fock basis states `sites_with_particle_superposition`
- superposition_coefficients are the coefficients of the superposition (assumed normalized)
"""
# TODO: fix, now broken 
function overlapwithfock(msum::MajoranaSum, sites_with_particle_superposition::Vector{FockState}, superposition_coefficients::Vector{<:Union{Real,Complex}})
    @error "Currently non supported"
    # check normalization
    @assert sum(abs2, superposition_coefficients) ≈ 1. "Superposition coefficients must be normalized."
    res = 0.
    unpaired_mask = create_unpaired_mask(nfermions(msum))

    for (ms, coeff) in zip(terms(msum), coefficients(msum))
        for (sites_with_particle, superposition_coefficient) in zip(sites_with_particle_superposition, superposition_coefficients)
            res += coeff * abs(superposition_coefficient)^2 * overlapwithfock(ms, unpaired_mask, sites_with_particle)
        end

        for k1 = 1:length(sites_with_particle_superposition)
            for k2 = k1+1:length(sites_with_particle_superposition)
                fock1 = sites_with_particle_superposition[k1]
                fock2 = sites_with_particle_superposition[k2]
                superposition_coeff1 = superposition_coefficients[k1]
                superposition_coeff2 = superposition_coefficients[k2]
                res += 2. * real(coeff * conj(superposition_coeff1) * superposition_coeff2 * overlapwithfock(ms, fock1, fock2, nfermions(msum)))
            end
        end
    end
    return res
end
