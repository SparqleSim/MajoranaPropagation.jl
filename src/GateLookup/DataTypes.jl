
mutable struct BaseCT{TT<:Integer}
    expression_pref::ComplexF64
    cos_theta_pref::Vector{Float64}
    sin_theta_pref::Vector{Float64}
    cumulative_string::TT
    function BaseCT(nfermions::Int)
        TT = getinttype(nfermions)
        return new{TT}(1., Float64[], Float64[], TT(0))
    end
    function BaseCT(expression_pref, cos_theta_pref, sin_theta_pref, cumulative_string::TT) where {TT<:Integer}
        return new{TT}(expression_pref, cos_theta_pref, sin_theta_pref, cumulative_string)
    end
end

mutable struct LookupCT{TT<:Integer}
    terms::Vector{BaseCT{TT}}
    function LookupCT(nfermions::Int)
        TT = getinttype(nfermions)
        return new{TT}([BaseCT(nfermions)])
    end
    function LookupCT(terms::Vector{BaseCT{TT}}) where {TT<:Integer}
        return new{TT}(terms)
    end
end

function LookupCT(::Type{TT}) where {TT<:Integer}
    return LookupCT(BaseCT{TT}[])
end

function Base.:+(coeff1::LookupCT{TT}, coeff2::LookupCT{TT}) where {TT<:Integer}
    return LookupCT(vcat(coeff1.terms, coeff2.terms))
end

