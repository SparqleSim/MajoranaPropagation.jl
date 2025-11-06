function order_sites(site_indices)
    sorted_indices = sort(site_indices)
    if sorted_indices != site_indices
        println("Warning: indices were not passed in ascending order")
        return sorted_indices
    end
    return site_indices
end

# TODO: list all supported symbols in the docstring below
# higher-level constructor for when passing symbols
# they wrap Symbol into Val for dispatch
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
function MajoranaSum(nfermions::Int, symbs::Vector{Symbol}, sites)
    obs = MajoranaSum(nfermions, Val(symbs[1]), sites[1])
    for i = 2:length(symbs)
        obs2 = MajoranaSum(nfermions, Val(symbs[i]), sites[i])
        obs = obs * obs2
    end
    return obs
end

function MajoranaSum(nfermions::Integer, symb::Symbol, sites)
    return MajoranaSum(nfermions, Val(symb), sites)
end

# Lower-level constructors for specific operators

# number operator
function MajoranaSum(nfermions::Integer, ::Val{:n}, site::Integer)
    TT = getinttype(nfermions)
    term1 = _bitonesat(TT, (2 * site - 1, 2 * site))
    term2 = TT(0)
    obs = MajoranaSum{TT,Float64}(nfermions, Dict(term1 => 0.5, term2 => 0.5))
    return obs
end

# hopping operator
function MajoranaSum(nfermions::Integer, ::Val{:h}, sites)
    TT = getinttype(nfermions)
    site1, site2 = order_sites(collect(sites))
    term1 = _bitonesat(TT, (2 * site1 - 1, 2 * site2))
    term2 = _bitonesat(TT, (2 * site1, 2 * site2 - 1))
    obs = MajoranaSum{TT,Float64}(nfermions, Dict(term1 => 0.5, term2 => -0.5))
    return obs
end

# number-number operator
function MajoranaSum(nfermions::Integer, ::Val{:nn}, sites)
    TT = getinttype(nfermions)
    site1, site2 = order_sites(collect(sites))
    term1 = _bitonesat(TT, (2 * site1 - 1, 2 * site1))
    term2 = _bitonesat(TT, (2 * site2 - 1, 2 * site2))
    term3 = _bitonesat(TT, (2 * site1 - 1, 2 * site1, 2 * site2 - 1, 2 * site2))
    term4 = TT(0)
    obs = MajoranaSum{TT,Float64}(
        nfermions,
        Dict(term1 => 0.25, term2 => 0.25, term3 => -0.25, term4 => 0.25)
    )
    return obs
end

# undefined
function MajoranaSum(nfermions::Integer, ::Val{symb}, sites) where {symb}
    error("Operator symbol :$symb not recognized.")
end


