
mutable struct FermionicRotationLookup{TT<:Integer} <: ParametrizedGate
    lookup_table::Dict{TT,Dict{TT,LookupCT{TT}}}
    n_sites::Integer
    is_spinful::Bool
    site_inds::Vector{Int}
end

function FermionicRotationLookup(symbol::Symbol, sites_acted_on::Integer, is_spinful::Bool, site_inds::Vector{Int}, theta::Float64)
    lookup_table = create_lookup_table(symbol, sites_acted_on, is_spinful)
    lookup_table = instantiate_lookup_table(lookup_table, theta)
    return FermionicRotationLookup(lookup_table, sites_acted_on, is_spinful, site_inds)
end



function create_lookup_table(symbol::Symbol, sites_acted_on, is_spinful)
    return create_lookup_table([symbol], sites_acted_on, is_spinful)
end

function create_lookup_table(symbols::Vector{Symbol}, sites_acted_on, is_spinful)
    gates = [FermionicRotation(symbol, collect(1:sites_acted_on)) for symbol in symbols]

    n_fermions = is_spinful ? 2 * sites_acted_on : sites_acted_on
    n_modes = 2 * n_fermions

    strings_to_check = MajoranaSum[]

    for k in 0:2:n_modes
        for combo in combinations(1:n_modes, k)
            msum = MajoranaSum(LookupCT, sites_acted_on, is_spinful)
            ms = MajoranaString(n_fermions, combo)
            add!(msum, ms.gammas, LookupCT(n_fermions))
            push!(strings_to_check, msum)
        end
    end

    @show strings_to_check
    @show length(strings_to_check)

    TT = getinttype(nfermions(strings_to_check[1]))

    lookup_table = Dict{TT,Dict{TT,LookupCT{TT}}}()

    for obs in strings_to_check
        println("------------------------------------------------")
        res = propagate(gates, obs, ones(length(gates)))
        @show obs, res 
        if length(res) > 1
            res_dict = Dict{TT,LookupCT{TT}}()
            for (ms, coeff) in res 
                res_dict[ms] = coeff 
            end
            lookup_table[collect(PropagationBase.terms(obs))[1]] = res_dict
        end
    end
    return lookup_table
end

function _all_prefactor_equal(lookup_table::Dict{TT,Dict{TT,LookupCT{TT}}}) where {TT<:Integer}
    return true 
end

function instantiate_coeff(coeff::BaseCT{TT}, cos_val::Float64, sin_val::Float64, cos2_val::Float64, sin2_val::Float64) where {TT<:Integer}
    new_coeff = BaseCT(coeff.expression_pref, [], [], coeff.cumulative_string)
    n_cos = length(coeff.cos_theta_pref)
    if n_cos == 0
    elseif n_cos == 1
        new_coeff.expression_pref *= cos_val 
    elseif n_cos == 2
        new_coeff.expression_pref *= cos2_val 
    else 
        new_coeff.expression_pref *= cos_val ^ n_cos
    end

    n_sin = length(coeff.sin_theta_pref)
    if n_sin == 0
    elseif n_sin == 1
        new_coeff.expression_pref *= sin_val 
    elseif n_sin == 2
        new_coeff.expression_pref *= sin2_val 
    else 
        new_coeff.expression_pref *= sin_val ^ n_sin
    end
    return new_coeff
end

function instantiate_coeff(coeff::LookupCT{TT}, cos_val::Float64, sin_val::Float64, cos2_val::Float64, sin2_val::Float64) where {TT<:Integer}
    new_coeff = LookupCT(TT)
    for k in eachindex(coeff.terms)
        push!(new_coeff.terms, instantiate_coeff(coeff.terms[k], cos_val, sin_val, cos2_val, sin2_val))
    end
    return new_coeff
end

function instantiate_lookup_table(lookup_table::Dict{TT,Dict{TT,LookupCT{TT}}}, theta) where {TT<:Integer}
    if _all_prefactor_equal(lookup_table)
        cos_val = cos(0.5 * theta)
        sin_val = sin(0.5 * theta)
        cos2_val = cos_val ^ 2
        sin2_val = sin_val ^ 2

        instantiated_table = Dict{TT,Dict{TT,LookupCT{TT}}}()

        for (k, v) in lookup_table
            v_inst = Dict{TT,LookupCT{TT}}()
            for (ms, coeff) in v 
                v_inst[ms] = instantiate_coeff(coeff, cos_val, sin_val, cos2_val, sin2_val)
            end
            instantiated_table[k] = v_inst
        end
    end
    return instantiated_table
end