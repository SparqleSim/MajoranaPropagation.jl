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
    MajoranaSum(n_sites::Integer, symb::Symbol, sites::Integer)
    MajoranaSum(n_sites::Integer, symb::Symbol, sites)

Returns a `MajoranaSum` corresponding to the observable defined by the symbol `symb` acting on one `site` or multiple `sites`.
Multiple sites can be passed as a vector or any iterable collection of integers.
The supported symbols are:
- Spinless operators 
    - `:n`: number operator on one given site
    - `:hop`: hopping operator between two given sites
    - `:nn`: number-number operator between the two given sites

- Spinful operators:
    - `:nup`: number operator for spin-up fermion on the given site
    - `:ndn`: number operator for spin-down fermion on the given site
    - `:nupndn`: number operator for undetermined spin fermion on the given site
    - `:hopup`: hopping operator for spin-up fermion between the two given sites
    - `:hopdn`: hopping operator for spin-down fermion between the two given sites
    - `:hole`: hole operator on the given site
"""

function MajoranaSum(n_sites::Integer, symb::Symbol, sites)
    return MajoranaSum(n_sites, Val(symb), sites)
end

# Lower-level constructors for specific operators

# Spinless operators

# number operator
function MajoranaSum(nfermions::Integer, ::Val{:n}, site)
    TT = getinttype(nfermions)
    is_spinful = false
    site = _tonum(site)
    term1 = _bitonesat(TT, (2 * site - 1, 2 * site))
    term2 = TT(0)
    obs = MajoranaSum{TT,Float64}(nfermions, is_spinful, Dict(term1 => 0.5, term2 => 0.5))
    return obs
end

# hopping operator
function MajoranaSum(nfermions::Integer, ::Val{:hop}, sites)
    TT = getinttype(nfermions)
    is_spinful = false
    sites = _tovec(sites)
    @assert length(sites) == 2 "Hopping operator requires exactly two site indices."
    site1, site2 = order_sites(sites)
    term1 = _bitonesat(TT, (2 * site1 - 1, 2 * site2))
    term2 = _bitonesat(TT, (2 * site1, 2 * site2 - 1))
    obs = MajoranaSum{TT,Float64}(nfermions, is_spinful, Dict(term1 => 0.5, term2 => -0.5))
    return obs
end

# number-number operator
function MajoranaSum(nfermions::Integer, ::Val{:nn}, sites)
    TT = getinttype(nfermions)
    is_spinful = false
    sites = _tovec(sites)
    @assert length(sites) == 2 "Number-number operator requires exactly two site indices."
    site1, site2 = order_sites(sites)
    term1 = _bitonesat(TT, (2 * site1 - 1, 2 * site1))
    term2 = _bitonesat(TT, (2 * site2 - 1, 2 * site2))
    term3 = _bitonesat(TT, (2 * site1 - 1, 2 * site1, 2 * site2 - 1, 2 * site2))
    term4 = TT(0)
    obs = MajoranaSum{TT,Float64}(
        nfermions,
        is_spinful,
        Dict(term1 => 0.25, term2 => 0.25, term3 => -0.25, term4 => 0.25)
    )
    return obs
end


# Spinful operators
function MajoranaSum(n_sites::Integer, ::Val{:nup}, site)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site = _tonum(site)
    term1 = _bitonesat(TT, (4 * site - 3, 4 * site - 2))
    term2 = TT(0)
    obs = MajoranaSum{TT,Float64}(n_sites, is_spinful, Dict(term1 => 0.5, term2 => 0.5))
    return obs
end

function MajoranaSum(n_sites::Integer, ::Val{:ndn}, site)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site = _tonum(site)
    term1 = _bitonesat(TT, (4 * site - 1, 4 * site))
    term2 = TT(0)
    obs = MajoranaSum{TT,Float64}(n_sites, is_spinful, Dict(term1 => 0.5, term2 => 0.5))
    return obs
end

function MajoranaSum(n_sites::Integer, ::Val{:hopup}, sites)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site1, site2 = order_sites(_tovec(sites))
    term1 = _bitonesat(TT, (4 * site1 - 3, 4 * site2 - 2))
    term2 = _bitonesat(TT, (4 * site1 - 2, 4 * site2 - 3))
    obs = MajoranaSum{TT,Float64}(n_sites, is_spinful, Dict(term1 => 0.5, term2 => -0.5))
    return obs
end

function MajoranaSum(n_sites::Integer, ::Val{:hopdn}, sites)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site1, site2 = order_sites(_tovec(sites))
    term1 = _bitonesat(TT, (4 * site1 - 1, 4 * site2))
    term2 = _bitonesat(TT, (4 * site1, 4 * site2 - 1))
    obs = MajoranaSum{TT,Float64}(n_sites, is_spinful, Dict(term1 => 0.5, term2 => -0.5))
    return obs
end

function MajoranaSum(n_sites::Integer, ::Val{:hole}, site)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site = _tonum(site)
    term1 = _bitonesat(TT, (4 * site - 3, 4 * site - 2))
    term2 = _bitonesat(TT, (4 * site - 1, 4 * site))
    term3 = _bitonesat(TT, (4 * site - 3, 4 * site - 2, 4 * site - 1, 4 * site))
    term4 = TT(0)
    obs = MajoranaSum{TT,Float64}(
        n_sites,
        is_spinful,
        Dict(term1 => -0.25, term2 => -0.25, term3 => -0.25, term4 => 0.25)
    )
    return obs
end

function MajoranaSum(n_sites::Integer, ::Val{:nupndn}, site)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site = _tonum(site)
    term1 = _bitonesat(TT, (4 * site - 3, 4 * site - 2))
    term2 = _bitonesat(TT, (4 * site - 1, 4 * site))
    term3 = _bitonesat(TT, (4 * site - 3, 4 * site - 2, 4 * site - 1, 4 * site))
    term4 = TT(0)
    obs = MajoranaSum{TT,Float64}(
        n_sites,
        is_spinful,
        Dict(term1 => 0.25, term2 => 0.25, term3 => -0.25, term4 => 0.25)
    )
    return obs
end

# undefined
function MajoranaSum(nfermions::Integer, ::Val{symb}, sites) where {symb}
    error("Operator symbol :$symb not recognized.")
end


_tovec(x) = collect(x)
_tovec(x::Number) = [x]
_tonum(x::Vector) = only(x)
_tonum(x::Number) = x