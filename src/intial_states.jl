function fockstate(n_sites::Integer, symb::Symbol, is_spinful::Bool; hole_positions=nothing, kwargs...)
    return fockstate(n_sites, Val(symb), is_spinful; hole_positions=hole_positions, kwargs...)
end

function fockstate(n_sites::Integer, ::Val{:checkerboard}, is_spinful::Bool; hole_positions=nothing, nx=-1, kwargs...)
    holes = isnothing(hole_positions) ? [] : hole_positions
    holes = isa(holes, Integer) ? [holes] : holes
    nx = nx == -1 ? n_sites : nx
    if is_spinful
        create_up_part_at = []
        create_down_part_at = []

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
        return fockstate(create_up_part_at, create_down_part_at)
    else
        @error "Checkerboard initial state is not supported for spinless fermions."
    end
    @error "Invalid symbol for fockstate constructor."
end

# printing for fock states
#=function Base.show(io::IO, fock_state::fockstate)
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
        occupied_sites = []
        for site in 1:fock_state.n_sites
            bit = (occupied_slots >> (site - 1)) & 1
            if bit == 1
                push!(occupied_sites, site)
            end
        end
        print(io, "occupied sites: ", occupied_sites)
    end
end=#
