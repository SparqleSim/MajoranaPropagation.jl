function order_sites(site_indices)
    sorted_indices = sort(site_indices)
    if sorted_indices != site_indices
        println("Warning: indices were not passed in ascending order")
        return sorted_indices
    end
    return site_indices
end

function ms_and_pref_from_symbol(nfermions::Int, symbs::Vector{Symbol}, sites, index::Int; pop_id = true)
    TT = getinttype(nfermions)
    symb = symbs[index]
    if symb == :n
        site = sites[index]
        obs = MajoranaSum(MajoranaString(nfermions, [2 * site - 1, 2 * site]), 0.5)
        if pop_id
            return obs, 0.5
        end 
        sum_add!(obs, TT(0), 0.5)
        return obs
    elseif symb == :h
        s1, s2 = order_sites(sites[index])
        obs = MajoranaSum(MajoranaString(nfermions, [2 * s1 - 1, 2 * s2]), 0.5)
        sum_add!(obs, MajoranaString(nfermions, [2 * s1, 2 * s2 - 1]), -0.5)
        if pop_id
            return obs, 0.0
        end
        return obs
    elseif symb == :nn 
        s1, s2 = order_sites(sites[index])
        obs = MajoranaSum(MajoranaString(nfermions, [2 * s1 - 1, 2 * s1]), 0.25)
        sum_add!(obs, MajoranaString(nfermions, [2 * s2 - 1, 2 * s2]), 0.25)
        sum_add!(obs, MajoranaString(nfermions, [2 * s1 - 1, 2 * s1, 2 * s2 - 1, 2 * s2]), -0.25)
        if pop_id
            return obs, 0.25
        end
        sum_add!(obs, TT(0), 0.25)
        return obs
    else
        throw(ArgumentError("Unknown symbol $symb for MajoranaSum"))
    end

end

function MajoranaSum(nfermions::Int, symbs::Vector{Symbol}, sites; pop_id = false)
    obs = ms_and_pref_from_symbol(nfermions, symbs, sites, 1; pop_id = false)
    for i = 2:length(symbs)
        obs2= ms_and_pref_from_symbol(nfermions, symbs, sites, i; pop_id = false)
        obs = obs * obs2 
    end
    if pop_id 
        TT = getinttype(nfermions)
        id_pref = pop!(obs, TT(0))
        return obs, id_pref
    end
    return obs
end

function MajoranaSum(nfermions::Int, symb::Symbol, sites; pop_id = false)
    return MajoranaSum(nfermions, [symb], [sites]; pop_id = pop_id)
end

