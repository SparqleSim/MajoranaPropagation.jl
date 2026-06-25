# ========================================================================== #
#  Building the Majorana transfer map for a fermionic gate.
#
#  The map is first built on a small "canonical" system whose sites are the gate's
#  distinct sites placed at 1, 2, ..., k. Each possible restricted Majorana string
#  on the gate's modes is propagated through the (trusted) `FermionicRotation`
#  implementation, giving its image under the gate's conjugate action. From the
#  image we extract, per output term, the cumulative gate string and a prefactor.
#
#  The map is then "shifted" onto the gate's actual modes (`site_inds`) so that
#  the cumulative strings live in the full system's bit layout and can be
#  multiplied against full observables with `ms_mult`.
# ========================================================================== #

"""
    _modes_per_site(is_spinful::Bool)

Number of Majorana modes per lattice site: 2 for spinless fermions (γ, γ'),
4 for spinful fermions (γ↑, γ'↑, γ↓, γ'↓).
"""
_modes_per_site(is_spinful::Bool) = is_spinful ? 4 : 2

"""
    _gate_modes(site_inds, is_spinful::Bool)

Sorted vector of (0-based) Majorana mode bit positions occupied by the gate.
Canonical mode `j` (1-based, into this vector) of the lookup table corresponds to
the actual bit position `_gate_modes(...)[j]`.
"""
function _gate_modes(site_inds, is_spinful::Bool)
    M = _modes_per_site(is_spinful)
    modes = Int[]
    for site in sort(collect(site_inds))
        for b in 0:(M-1)
            push!(modes, (site - 1) * M + b)
        end
    end
    return modes
end

"""
    _compress(majstring, gate_modes)

Extract the bits of `majstring` at positions `gate_modes` and pack them into the
low `length(gate_modes)` bits, giving the 0-based column index into a transfer map.
"""
@inline function _compress(majstring::Integer, gate_modes::Vector{Int})
    idx = 0
    @inbounds for j in eachindex(gate_modes)
        idx |= Int((majstring >> gate_modes[j]) & 1) << (j - 1)
    end
    return idx
end

"""
    _expand(canonical_string, gate_modes, ::Type{TT})

Inverse of [`_compress`](@ref): place the low bits of `canonical_string` onto the
actual mode positions `gate_modes`, returning a string of type `TT`.
"""
function _expand(canonical_string::Integer, gate_modes::Vector{Int}, ::Type{TT}) where {TT<:Integer}
    result = zero(TT)
    @inbounds for j in eachindex(gate_modes)
        if (canonical_string >> (j - 1)) & 1 == 1
            result |= TT(1) << gate_modes[j]
        end
    end
    return result
end

"""
    _build_raw_columns(symbols, n_canonical_sites, canonical_sites, is_spinful, theta)

Propagate every restricted Majorana string on the gate's modes through the
canonical fermionic gate(s) and return `(columns, n_fermions_canonical)`.

`n_canonical_sites` is the number of distinct sites of the canonical system (`1:n_canonical_sites`).
`canonical_sites` are the rank-mapped sites the canonical gate acts on (values in
`1:n_canonical_sites`, possibly repeated); passing them in the original site order preserves
direction-sensitive gates.

`columns[s+1]` is the image of the (canonical) input string `s`: a list of
`(output_string, coeff)` tuples with `output_string` and `coeff` in the canonical frame.
"""
function _build_raw_columns(symbols::Vector{Symbol}, n_canonical_sites::Integer, canonical_sites::Vector{Int}, is_spinful::Bool, theta::Real)
    M = _modes_per_site(is_spinful)
    n_modes = n_canonical_sites * M
    n_fermions = is_spinful ? 2 * n_canonical_sites : n_canonical_sites
    TT = getinttype(n_fermions)

    gates = [FermionicRotation(symbol, canonical_sites) for symbol in symbols]
    thetas = fill(Float64(theta), length(gates))

    columns = Vector{Vector{Tuple{TT,ComplexF64}}}(undef, 2^n_modes)

    for s in 0:(2^n_modes - 1)
        # single-term observable: the restricted Majorana string `s` with coefficient 1
        obs = MajoranaSum(Float64, n_canonical_sites, is_spinful)
        add!(obs, TT(s), 1.0)

        # propagate through the gate(s); keep all terms (no truncation)
        res = propagate(gates, obs, thetas; min_abs_coeff=-1.0)

        column = Tuple{TT,Float64}[]
        for (out_string, coeff) in res
            if abs(coeff) < 1e-12
                continue
            end
            push!(column, (TT(out_string), coeff))
        end
        # an empty column means the gate acts as the identity on this string
        if isempty(column)
            push!(column, (TT(s), 1.0))
        end
        columns[s+1] = column
    end

    return columns, n_fermions
end

"""
    _finalize_columns(raw_columns, ::Type{TTin}, transform, n_fermions, ::Type{TTout})

Convert the raw `(output_string, coeff)` columns into transfer-map columns of
`(cumulative_string, ctilde)` tuples expressed in the frame defined by `transform`.

For an input string `s` and output `(s_out, a)`, the cumulative gate string is
`C = transform(s) ⊻ transform(s_out)` and `ctilde = a / μ`, where
`μ = ms_mult(C, transform(s))` is the sign of multiplying the cumulative string onto
the (restricted) input in this frame. At application time the actual coefficient of the
output string is `real(η * ctilde)`, with `η` the sign of `ms_mult(C, full_observable)`;
the `μ`/`η` split isolates the full-string-dependent part into `η`.
"""
function _finalize_columns(raw_columns::Vector{Vector{Tuple{TTin,ComplexF64}}}, ::Type{TTin},
    transform, n_fermions::Integer, ::Type{TTout}) where {TTin<:Integer,TTout<:Integer}

    columns = Vector{Vector{Tuple{TTout,ComplexF64}}}(undef, length(raw_columns))
    for i in eachindex(raw_columns)
        s_in = transform(TTin(i - 1))::TTout
        column = Tuple{TTout,ComplexF64}[]
        for (s_out, a) in raw_columns[i]
            s_out_t = transform(s_out)::TTout
            cumulative = s_in ⊻ s_out_t
            mu, _ = ms_mult(cumulative, s_in, n_fermions)
            push!(column, (cumulative, a / mu))
        end
        columns[i] = column
    end
    return columns
end

# ========================================================================== #
#  Angle-free (surrogate-backed) table building.
#
#  Identical in structure to `_build_raw_columns` / `_finalize_columns`, but each
#  basis string is propagated through the *surrogate* (see `src/surrogate.jl`) so
#  that the propagated coefficient `a` is kept as a node graph in the rotation angle
#  rather than a number. The angle-independent structure (cumulative strings and the
#  frame prefactor `1/μ`) is computed once at build time; `evaluate` later walks the
#  node graphs to obtain `a` and forms `ctilde = a / μ` for a concrete angle.
# ========================================================================== #

"""
    _build_symbolic_columns(symbols, n_canonical_sites, canonical_sites, is_spinful)

Surrogate analog of [`_build_raw_columns`](@ref). Propagates every restricted basis Majorana
string through the canonical gate(s) using the surrogate and returns
`(columns, n_fermions_canonical, nparams)`, where `columns[s+1]` is a list of
`(output_string, path)` tuples (`path::MajoranaNodePathProperties` carrying the symbolic angle
dependence) and `nparams = length(symbols)` is the number of angles `evaluate` must supply.
All structurally-produced terms are kept (no coefficient filtering): a term that vanishes at one
angle may be nonzero at another.
"""
function _build_symbolic_columns(symbols::Vector{Symbol}, n_canonical_sites::Integer, canonical_sites::Vector{Int}, is_spinful::Bool)
    M = _modes_per_site(is_spinful)
    n_modes = n_canonical_sites * M
    n_fermions = is_spinful ? 2 * n_canonical_sites : n_canonical_sites
    TT = getinttype(n_fermions)

    gates = [FermionicRotation(symbol, canonical_sites) for symbol in symbols]
    nparams = length(gates)

    columns = Vector{Vector{Tuple{TT,MajoranaNodePathProperties}}}(undef, 2^n_modes)

    for s in 0:(2^n_modes - 1)
        # single-term observable: the restricted Majorana string `s` with coefficient 1
        obs = MajoranaSum(Float64, n_canonical_sites, is_spinful)
        add!(obs, TT(s), 1.0)
        obs_surr = wrapcoefficients(obs, MajoranaNodePathProperties)

        # propagate through the surrogate (no thetas; builds the symbolic graph once)
        res = propagate(gates, obs_surr)

        column = Tuple{TT,MajoranaNodePathProperties}[]
        for (out_string, path) in res
            push!(column, (TT(out_string), path))
        end
        # a gate acting as the identity leaves the input string with a constant coefficient;
        # this fallback is defensive (the surrogate never truncates, so `res` is non-empty)
        if isempty(column)
            push!(column, (TT(s), first(coefficients(obs_surr))))
        end
        columns[s+1] = column
    end

    return columns, n_fermions, nparams
end

"""
    _finalize_symbolic_columns(raw_columns, ::Type{TTin}, transform, n_fermions, ::Type{TTout})

Surrogate analog of [`_finalize_columns`](@ref). Computes the angle-independent
`(cumulative_string, SurrogateCoeff(path, 1/μ))` entries in the frame defined by `transform`,
leaving the symbolic coefficient `path` untouched (the node graphs are shared between the
canonical and shifted maps, so they are evaluated only once).
"""
function _finalize_symbolic_columns(raw_columns::Vector{Vector{Tuple{TTin,MajoranaNodePathProperties}}}, ::Type{TTin},
    transform, n_fermions::Integer, ::Type{TTout}) where {TTin<:Integer,TTout<:Integer}

    columns = Vector{Vector{Tuple{TTout,SurrogateCoeff}}}(undef, length(raw_columns))
    for i in eachindex(raw_columns)
        s_in = transform(TTin(i - 1))::TTout
        column = Tuple{TTout,SurrogateCoeff}[]
        for (s_out, path) in raw_columns[i]
            s_out_t = transform(s_out)::TTout
            cumulative = s_in ⊻ s_out_t
            mu, _ = ms_mult(cumulative, s_in, n_fermions)
            push!(column, (cumulative, SurrogateCoeff(path, ComplexF64(1 / mu))))
        end
        columns[i] = column
    end
    return columns
end
