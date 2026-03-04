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
