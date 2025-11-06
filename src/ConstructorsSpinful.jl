# just like in Constructors.jl, but for spinful operators
# we have the high level there so only define Val dispatch here

# spinful operators on n sites have 2*n fermions

function MajoranaSum(spinful_sites::Integer, ::Val{:nu}, site::Integer)
    nfermions = 2 * spinful_sites
    TT = getinttype(nfermions)
    term1 = _bitonesat(TT, (4 * site - 3, 4 * site - 2))
    term2 = TT(0)
    obs = MajoranaSum{TT,Float64}(nfermions, Dict(term1 => 0.5, term2 => 0.5))
    return obs
end

function MajoranaSum(spinful_sites::Integer, ::Val{:nd}, site::Integer)
    nfermions = 2 * spinful_sites
    TT = getinttype(nfermions)
    term1 = _bitonesat(TT, (4 * site - 1, 4 * site))
    term2 = TT(0)
    obs = MajoranaSum{TT,Float64}(nfermions, Dict(term1 => 0.5, term2 => 0.5))
    return obs
end

function MajoranaSum(spinful_sites::Integer, ::Val{:hu}, sites)
    nfermions = 2 * spinful_sites
    TT = getinttype(nfermions)
    site1, site2 = order_sites(collect(sites))
    term1 = _bitonesat(TT, (4 * site1 - 3, 4 * site2 - 2))
    term2 = _bitonesat(TT, (4 * site1 - 2, 4 * site2 - 3))
    obs = MajoranaSum{TT,Float64}(nfermions, Dict(term1 => 0.5, term2 => -0.5))
    return obs
end

function MajoranaSum(spinful_sites::Integer, ::Val{:hd}, sites)
    nfermions = 2 * spinful_sites
    TT = getinttype(nfermions)
    site1, site2 = order_sites(collect(sites))
    term1 = _bitonesat(TT, (4 * site1 - 1, 4 * site2))
    term2 = _bitonesat(TT, (4 * site1, 4 * site2 - 1))
    obs = MajoranaSum{TT,Float64}(nfermions, Dict(term1 => 0.5, term2 => -0.5))
    return obs
end

function MajoranaSum(spinful_sites::Integer, ::Val{:hole}, site::Integer)
    nfermions = 2 * spinful_sites
    TT = getinttype(nfermions)
    term1 = _bitonesat(TT, (4 * site - 3, 4 * site - 2))
    term2 = _bitonesat(TT, (4 * site - 1, 4 * site))
    term3 = _bitonesat(TT, (4 * site - 3, 4 * site - 2, 4 * site - 1, 4 * site))
    term4 = TT(0)
    obs = MajoranaSum{TT,Float64}(
        nfermions,
        Dict(term1 => -0.25, term2 => -0.25, term3 => -0.25, term4 => 0.25)
    )
    return obs
end

function MajoranaSum(spinful_sites::Integer, ::Val{:nund}, site::Integer)
    nfermions = 2 * spinful_sites
    TT = getinttype(nfermions)
    term1 = _bitonesat(TT, (4 * site - 3, 4 * site - 2))
    term2 = _bitonesat(TT, (4 * site - 1, 4 * site))
    term3 = _bitonesat(TT, (4 * site - 3, 4 * site - 2, 4 * site - 1, 4 * site))
    term4 = TT(0)
    obs = MajoranaSum{TT,Float64}(
        nfermions,
        Dict(term1 => 0.25, term2 => 0.25, term3 => -0.25, term4 => 0.25)
    )
    return obs
end
