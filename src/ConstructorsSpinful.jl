function spinful_ms_and_pref_from_symbol(n_spinful_sites::Int, symbs::Vector{Symbol}, sites, index::Int; pop_id = true)
    TT = getinttype(2 * n_spinful_sites)
    symb = symbs[index]
    if symb == :nu
        site = sites[index]
        obs = MajoranaSum(MajoranaString(2 * n_spinful_sites, [4 * site - 3, 4 * site - 2]), 0.5)

        if pop_id
            return obs, 0.5
        end

        sum_add!(obs, TT(0), 0.5)
        return obs
    elseif symb == :nd
        site = sites[index]
        obs =  MajoranaSum(MajoranaString(2 * n_spinful_sites, [4 * site - 1, 4 * site]), 0.5)

        if pop_id
            return obs, 0.5
        end

        sum_add!(obs, TT(0), 0.5)
        return obs
    elseif symb == :hu
        s1, s2 = order_sites(sites[index])
        obs = MajoranaSum(MajoranaString(2 * n_spinful_sites, [4 * s1 - 3, 4 * s2 - 2]), 0.5)
        sum_add!(obs, MajoranaString(2 * n_spinful_sites, [4 * s1 - 2, 4 * s2 - 3]), -0.5)
        if pop_id
            return obs, 0.0
        end
        return obs
    elseif symb == :hd
        s1, s2 = order_sites(sites[index])
        obs = MajoranaSum(MajoranaString(2 * n_spinful_sites, [4 * s1 - 1, 4 * s2]), 0.5)
        sum_add!(obs, MajoranaString(2 * n_spinful_sites, [4 * s1, 4 * s2 - 1]), -0.5)
        if pop_id
            return obs, 0.0
        end
        return obs
    elseif symb == :hole 
        site = sites[index]
        obs = MajoranaSum(MajoranaString(2 * n_spinful_sites, [4 * site - 3, 4 * site - 2]), -0.25)
        sum_add!(obs, MajoranaString(2 * n_spinful_sites, [4 * site - 1, 4 * site]), -0.25)
        sum_add!(obs, MajoranaString(2 * n_spinful_sites, [4 * site - 3, 4 * site - 2, 4 * site - 1, 4 * site]), -0.25)
        if pop_id
            return obs, 0.25
        end
        sum_add!(obs, TT(0), 0.25)
        return obs
    elseif symb == :nund 
        site = sites[index]
        obs = MajoranaSum(MajoranaString(2 * n_spinful_sites, [4 * site - 3, 4 * site - 2]), 0.25)
        sum_add!(obs, MajoranaString(2 * n_spinful_sites, [4 * site - 1, 4 * site]), 0.25)
        sum_add!(obs, MajoranaString(2 * n_spinful_sites, [4 * site - 3, 4 * site - 2, 4 * site - 1, 4 * site]), -0.25)
        if pop_id
            return obs, 0.25
        end
        sum_add!(obs, TT(0), 0.25)
        return obs
    else 
        throw(ArgumentError("Unknown symbol $symb for MajoranaSum"))
    end
end

function spinfulmajoranasum(n_spinful_sites::Int, symbs::Vector{Symbol}, sites; pop_id = false)
    @assert length(symbs) == length(sites)
    obs = spinful_ms_and_pref_from_symbol(n_spinful_sites, symbs, sites, 1; pop_id = false)
    for i = 2:length(symbs)
        obs2= spinful_ms_and_pref_from_symbol(n_spinful_sites, symbs, sites, i; pop_id = false)
        obs = obs * obs2
    end

    if pop_id 
        TT = getinttype(2 * n_spinful_sites)
        id_pref = pop!(obs, TT(0))
        return obs, id_pref
    end
    return obs
end

function spinfulmajoranasum(n_spinful_sites::Int, symb, sites; pop_id = false)
    return spinfulmajoranasum(n_spinful_sites, [symb], [sites]; pop_id = pop_id)
end