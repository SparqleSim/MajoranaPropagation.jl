function order_sites(site_indices)
    sorted_indices = sort(site_indices)
    if sorted_indices != site_indices
        @warn "Indices were not passed in ascending order, which is always assumed by the constructor"
        return sorted_indices
    end
    return site_indices
end

# higher-level constructor for when passing symbols
# they wrap Symbol into Val for dispatch
"""
    MajoranaSum(n_sites::Integer, symb::Symbol, sites)

Return a `MajoranaSum` corresponding to the observable defined by the symbol `symb` acting on one site or multiple `sites` (a single integer, or any iterable of integers).
For spinless symbols the first argument counts the fermions, for spinful symbols it counts the (two-fermion) sites.
Two-site operators are defined with the site indices sorted as ``i < j``, except for `:hopupdn` (see below).
The fermionic annihilation and creation operators on site ``i`` are denoted ``f_i, f_i^\\dagger``, with ``n_i = f_i^\\dagger f_i``.

Spinless operators:
- `:n`: number operator ``n_i`` on the given site
- `:hop`: hopping operator ``f_i^\\dagger f_j + f_j^\\dagger f_i`` between the two given sites
- `:nn`: density-density operator ``n_i n_j`` between the two given sites
- `:pair`: pairing operator ``f_i^\\dagger f_j^\\dagger + f_j f_i`` between the two given sites
- `:f`: annihilation operator ``f_i`` on the given site (complex coefficients)
- `:fdag`: creation operator ``f_i^\\dagger`` on the given site (complex coefficients)

Spinful operators:
- `:nup`, `:ndn`: number operator ``n_{i\\uparrow}`` / ``n_{i\\downarrow}`` on the given site
- `:nupndn`: double-occupancy operator ``n_{i\\uparrow} n_{i\\downarrow}`` on the given site
- `:hole`: hole operator ``(1 - n_{i\\uparrow})(1 - n_{i\\downarrow})`` on the given site
- `:hopup`, `:hopdn`: hopping operator ``f_{i\\uparrow}^\\dagger f_{j\\uparrow} + f_{j\\uparrow}^\\dagger f_{i\\uparrow}`` (resp. spin down) between the two given sites
- `:hop_on_site`: on-site spin-flip operator ``f_{i\\uparrow}^\\dagger f_{i\\downarrow} + f_{i\\downarrow}^\\dagger f_{i\\uparrow}`` on the given site
- `:hopupdn`: spin-mixing hopping ``f_{i\\uparrow}^\\dagger f_{j\\downarrow} + f_{j\\downarrow}^\\dagger f_{i\\uparrow}``, where the first given site carries the spin-up and the second the spin-down fermion (the two sites are not sorted)
- `:pairup`, `:pairdn`: pairing operator ``f_{i\\uparrow}^\\dagger f_{j\\uparrow}^\\dagger + f_{j\\uparrow} f_{i\\uparrow}`` (resp. spin down) between the two given sites
- `:fup`, `:fupdag`, `:fdn`, `:fdndag`: annihilation / creation operators ``f_{i\\uparrow}``, ``f_{i\\uparrow}^\\dagger``, ``f_{i\\downarrow}``, ``f_{i\\downarrow}^\\dagger`` on the given site (complex coefficients)
- `:Sz`: spin operator ``(n_{i\\uparrow} - n_{i\\downarrow})/2`` on the given site
"""
function MajoranaSum(n_sites::Integer, symb::Symbol, sites)
    return MajoranaSum(n_sites, Val(symb), sites)
end

# Lower-level constructors for specific operators

# ========== Spinless operators ==========

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

#pair operator
function MajoranaSum(nfermions::Integer, ::Val{:pair}, sites)
    TT = getinttype(nfermions)
    is_spinful = false
    site1, site2 = order_sites(_tovec(sites))
    term1 = _bitonesat(TT, (2 * site1 - 1, 2 * site2))
    term2 = _bitonesat(TT, (2 * site1, 2 * site2 - 1))
    obs = MajoranaSum{TT,Float64}(nfermions, is_spinful, Dict(term1 => -0.5, term2 => -0.5))
    return obs
end

# annhililation operator
function MajoranaSum(nfermions::Integer, ::Val{:f}, site)
    TT = getinttype(nfermions)
    is_spinful = false
    site = _tonum(site)
    term1 = _bitonesat(TT, (2 * site - 1))
    term2 = _bitonesat(TT, (2 * site))
    obs = MajoranaSum{TT,ComplexF64}(nfermions, is_spinful, Dict(term1 => 0.5 + 0.0im, term2 => 0.5im))
    return obs
end
# creation operator
function MajoranaSum(nfermions::Integer, ::Val{:fdag}, site)
    TT = getinttype(nfermions)
    is_spinful = false
    site = _tonum(site)
    term1 = _bitonesat(TT, (2 * site - 1))
    term2 = _bitonesat(TT, (2 * site))
    obs = MajoranaSum{TT,ComplexF64}(nfermions, is_spinful, Dict(term1 => 0.5 + 0.0im, term2 => -0.5im))
    return obs
end

# ========== Spinful operators ==========

# ========== Spinful operators ==========

# n_up operator
function MajoranaSum(n_sites::Integer, ::Val{:nup}, site)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site = _tonum(site)
    term1 = _bitonesat(TT, (4 * site - 3, 4 * site - 2))
    term2 = TT(0)
    obs = MajoranaSum{TT,Float64}(n_sites, is_spinful, Dict(term1 => 0.5, term2 => 0.5))
    return obs
end

# n_dn operator
function MajoranaSum(n_sites::Integer, ::Val{:ndn}, site)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site = _tonum(site)
    term1 = _bitonesat(TT, (4 * site - 1, 4 * site))
    term2 = TT(0)
    obs = MajoranaSum{TT,Float64}(n_sites, is_spinful, Dict(term1 => 0.5, term2 => 0.5))
    return obs
end

# up hopping operator
function MajoranaSum(n_sites::Integer, ::Val{:hopup}, sites)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site1, site2 = order_sites(_tovec(sites))
    term1 = _bitonesat(TT, (4 * site1 - 3, 4 * site2 - 2))
    term2 = _bitonesat(TT, (4 * site1 - 2, 4 * site2 - 3))
    obs = MajoranaSum{TT,Float64}(n_sites, is_spinful, Dict(term1 => 0.5, term2 => -0.5))
    return obs
end

# down hopping operator
function MajoranaSum(n_sites::Integer, ::Val{:hopdn}, sites)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site1, site2 = order_sites(_tovec(sites))
    term1 = _bitonesat(TT, (4 * site1 - 1, 4 * site2))
    term2 = _bitonesat(TT, (4 * site1, 4 * site2 - 1))
    obs = MajoranaSum{TT,Float64}(n_sites, is_spinful, Dict(term1 => 0.5, term2 => -0.5))
    return obs
end

# onsite hopping operator
function MajoranaSum(n_sites::Integer, ::Val{:hop_on_site}, site)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site = _tonum(site)
    term1 = _bitonesat(TT, (4 * site - 3, 4 * site))
    term2 = _bitonesat(TT, (4 * site - 2, 4 * site - 1))
    obs = MajoranaSum{TT,Float64}(n_sites, is_spinful, Dict(term1 => 0.5, term2 => -0.5))
    return obs
end

# up-down hopping operator
# up: site[1]
# down: site[2]
function MajoranaSum(n_sites::Integer, ::Val{:hopupdn}, sites)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site1, site2 = _tovec(sites) # here it's important to not order the sites, since they refer to different spins
    term1 = _bitonesat(TT, (4 * site1 - 3, 4 * site2))
    term2 = _bitonesat(TT, (4 * site1 - 2, 4 * site2 - 1))
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

function MajoranaSum(n_sites::Integer, ::Val{:pairup}, sites)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site1, site2 = order_sites(_tovec(sites))
    term1 = _bitonesat(TT, (4 * site1 - 3, 4 * site2 - 2))
    term2 = _bitonesat(TT, (4 * site1 - 2, 4 * site2 - 3))
    obs = MajoranaSum{TT,Float64}(n_sites, is_spinful, Dict(term1 => -0.5, term2 => -0.5))
    return obs
end

function MajoranaSum(n_sites::Integer, ::Val{:pairdn}, sites)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site1, site2 = order_sites(_tovec(sites))
    term1 = _bitonesat(TT, (4 * site1 - 1, 4 * site2))
    term2 = _bitonesat(TT, (4 * site1, 4 * site2 - 1))
    obs = MajoranaSum{TT,Float64}(n_sites, is_spinful, Dict(term1 => -0.5, term2 => -0.5))
    return obs
end

# up annihilation operator
function MajoranaSum(n_sites::Integer, ::Val{:fup}, site)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site = _tonum(site)
    term1 = _bitonesat(TT, (4 * site - 3))
    term2 = _bitonesat(TT, (4 * site - 2))
    obs = MajoranaSum{TT,ComplexF64}(n_sites, is_spinful, Dict(term1 => 0.5 + 0.0im, term2 => 0.5im))
    return obs
end
# up creation operator
function MajoranaSum(n_sites::Integer, ::Val{:fupdag}, site)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site = _tonum(site)
    term1 = _bitonesat(TT, (4 * site - 3))
    term2 = _bitonesat(TT, (4 * site - 2))
    obs = MajoranaSum{TT,ComplexF64}(n_sites, is_spinful, Dict(term1 => 0.5 + 0.0im, term2 => -0.5im))
    return obs
end
# down annihilation operator
function MajoranaSum(n_sites::Integer, ::Val{:fdn}, site)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site = _tonum(site)
    term1 = _bitonesat(TT, (4 * site - 1))
    term2 = _bitonesat(TT, (4 * site))
    obs = MajoranaSum{TT,ComplexF64}(n_sites, is_spinful, Dict(term1 => 0.5 + 0.0im, term2 => 0.5im))
    return obs
end
# down creation operator
function MajoranaSum(n_sites::Integer, ::Val{:fdndag}, site)
    TT = getinttype(2 * n_sites)
    is_spinful = true
    site = _tonum(site)
    term1 = _bitonesat(TT, (4 * site - 1))
    term2 = _bitonesat(TT, (4 * site))
    obs = MajoranaSum{TT,ComplexF64}(n_sites, is_spinful, Dict(term1 => 0.5 + 0.0im, term2 => -0.5im))
    return obs
end

#Sz operator 
function MajoranaSum(n_sites::Integer, ::Val{:Sz}, site)
    msum = 0.5 * MajoranaSum(n_sites, :nup, site) - 0.5 * MajoranaSum(n_sites, :ndn, site)
    delete!(msum, 0)
    return msum
end

# undefined
function MajoranaSum(nfermions::Integer, ::Val{symb}, sites) where {symb}
    error("Operator symbol :$symb not recognized.")
end


_tovec(x) = collect(x)
_tovec(x::Number) = [x]
_tonum(x::Vector) = only(x)
_tonum(x::Number) = x


"""
    flag_non_number_preserving(symb::Symbol)
Given a fermionic gate symbol, returns true if the decomposition of the gate into
Majorana rotations would lead to an unphysical non-number-preserving operation (e.g., hoppings), false otherwise.
This is then used to determine whether to truncate after each Majorana rotation when applying the gate.
"""
function flag_non_number_preserving(symb::Symbol)
    return flag_non_number_preserving(Val(symb))
end

# default val, don't flag
function flag_non_number_preserving(::Val{symb}) where {symb}
    return false
end
#flag hoppings
function flag_non_number_preserving(::Val{:hop})
    return true
end
function flag_non_number_preserving(::Val{:hopup})
    return true
end
function flag_non_number_preserving(::Val{:hopdn})
    return true
end
