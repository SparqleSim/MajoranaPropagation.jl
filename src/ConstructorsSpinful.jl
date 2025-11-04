function spinful_ms_and_pref_from_symbol(n_spinful_sites::Int, symbs::Vector{Symbol}, sites, index::Int)
    TT = getinttype(2 * n_spinful_sites)
    symb = symbs[index]
    if symb == :nu
        site = sites[index]
        obs = MajoranaSum(MajoranaString(2 * n_spinful_sites, [4 * site - 3, 4 * site - 2]), 0.5)
        sum_add!(obs, TT(0), 0.5)
        return obs
    elseif symb == :nd
        site = sites[index]
        obs =  MajoranaSum(MajoranaString(2 * n_spinful_sites, [4 * site - 1, 4 * site]), 0.5)
        sum_add!(obs, TT(0), 0.5)
        return obs
    elseif symb == :hu
        s1, s2 = order_sites(sites[index])
        obs = MajoranaSum(MajoranaString(2 * n_spinful_sites, [4 * s1 - 3, 4 * s2 - 2]), 0.5)
        sum_add!(obs, MajoranaString(2 * n_spinful_sites, [4 * s1 - 2, 4 * s2 - 3]), -0.5)
        return obs
    elseif symb == :hd
        s1, s2 = order_sites(sites[index])
        obs = MajoranaSum(MajoranaString(2 * n_spinful_sites, [4 * s1 - 1, 4 * s2]), 0.5)
        sum_add!(obs, MajoranaString(2 * n_spinful_sites, [4 * s1, 4 * s2 - 1]), -0.5)
        return obs
    elseif symb == :hole 
        site = sites[index]
        obs = MajoranaSum(MajoranaString(2 * n_spinful_sites, [4 * site - 3, 4 * site - 2]), -0.25)
        sum_add!(obs, MajoranaString(2 * n_spinful_sites, [4 * site - 1, 4 * site]), -0.25)
        sum_add!(obs, MajoranaString(2 * n_spinful_sites, [4 * site - 3, 4 * site - 2, 4 * site - 1, 4 * site]), -0.25)
        sum_add!(obs, TT(0), 0.25)
        return obs
    elseif symb == :nund 
        site = sites[index]
        obs = MajoranaSum(MajoranaString(2 * n_spinful_sites, [4 * site - 3, 4 * site - 2]), 0.25)
        sum_add!(obs, MajoranaString(2 * n_spinful_sites, [4 * site - 1, 4 * site]), 0.25)
        sum_add!(obs, MajoranaString(2 * n_spinful_sites, [4 * site - 3, 4 * site - 2, 4 * site - 1, 4 * site]), -0.25)
        sum_add!(obs, TT(0), 0.25)
        return obs
    else 
        throw(ArgumentError("Unknown symbol $symb for MajoranaSum"))
    end
end

""" 
    spinfulmajoranasum(n_spinful_sites::Int, symbs::Vector{Symbol}, sites)

Returns a `MajoranaSum` with with 2*n_spinful_sites fermionic sites, corresponding to the observable defined by the product
of the symbols in `symbs` acting on the sites in `sites`.
The supported symbols are:
- `:nu`: number operator for spin-up fermion on the given site
- `:nd`: number operator for spin-down fermion on the given site
- `:nund`: number operator for undetermined spin fermion on the given site
- `:hu`: hopping operator for spin-up fermion between the two given sites
- `:hd`: hopping operator for spin-down fermion between the two given sites
- `:hole`: hole operator on the given site
The sites should be given as lists. For example 
`sites = [[1, 2]]` for a symbol acting on site 1 and site 2 respectively.

They can be combined, eg. nu_1 * nu_2 corresponds to `symbs = [:nu, :nu]` and `sites = [1, 2]`.
"""
function spinfulmajoranasum(n_spinful_sites::Int, symbs::Vector{Symbol}, sites)
    @assert length(symbs) == length(sites)
    obs = spinful_ms_and_pref_from_symbol(n_spinful_sites, symbs, sites, 1)
    for i = 2:length(symbs)
        obs2= spinful_ms_and_pref_from_symbol(n_spinful_sites, symbs, sites, i)
        obs = obs * obs2
    end
    return obs
end


""" 
    spinfulmajoranasum(n_spinful_sites::Int, symbs::Vector{Symbol}, sites)

Returns a `MajoranaSum` with with 2*n_spinful_sites fermionic sites, corresponding to the observable defined by
the symbol in `symb` acting on the sites in `sites`.
"""
function spinfulmajoranasum(n_spinful_sites::Int, symb, sites)
    return spinfulmajoranasum(n_spinful_sites, [symb], [sites])
end