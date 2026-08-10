# Dense Jordan-Wigner reference implementation ("oracle") for small systems.
#
# Everything in here builds explicit 2^nf x 2^nf matrices, so it is only usable
# for small mode counts (nf <= ~10). It provides ground truth for overlaps,
# products, gates and imaginary-time evolution, independent of the bit-twiddling
# in src/.
#
# Conventions (verified against the package, see test_fock_dense.jl which
# cross-validates this oracle against the Yao-checked compare_jw.jl path):
#
# - Majorana mode index k = 1..2nf; fermionic site s carries gamma_s at
#   k = 2s-1 and gamma'_s at k = 2s (src/MajoranaDataTypes.jl). For spinful
#   systems site s maps to modes 2s-1 (up) and 2s (down), so the same oracle
#   works with nf = 2 * n_sites.
# - Jordan-Wigner: gamma_{2s-1} = Z^(s-1) X_s, gamma_{2s} = Z^(s-1) Y_s,
#   with vacuum |0...0> and n_s = (1 - Z_s)/2. Then a_s = (gamma_s + i gamma'_s)/2.
# - A stored (string, coeff) pair of a MajoranaSum represents the Hermitian
#   operator coeff * i^omega_L_mult(ms) * prod_{k ascending} gamma_k,
#   where omega_L_mult(ms) = w(w-1)/2 mod 2 for weight w.
#   (E.g. :n on site s is stored as {gamma_s gamma'_s => 0.5, I => 0.5},
#   representing 0.5 * i gamma gamma' + 0.5 * I = a'a.)
# - FockState phase convention: |f> is the *descending* product of creation
#   operators over occupied modes, i.e. eta_f = (-1)^(n(n-1)/2) relative to the
#   ascending-product basis vector, with n the particle number. This matters
#   only for off-diagonal matrix elements and superposition overlaps.

using LinearAlgebra

const _X = ComplexF64[0 1; 1 0]
const _Y = ComplexF64[0 -im; im 0]
const _Z = ComplexF64[1 0; 0 -1]
const _I2 = ComplexF64[1 0; 0 1]

"""
    dense_majoranas(nf)
Return the 2nf dense Majorana matrices gamma_1, ..., gamma_2nf via Jordan-Wigner.
"""
function dense_majoranas(nf::Int)
    gammas = Matrix{ComplexF64}[]
    for k in 1:2nf
        s = (k + 1) ÷ 2
        ops = [j < s ? _Z : (j == s ? (isodd(k) ? _X : _Y) : _I2) for j in 1:nf]
        push!(gammas, reduce(kron, ops))
    end
    return gammas
end

"""
    dense_string(ms_int, nf, gammas)
Dense matrix of the stored Majorana string `ms_int` (an integer bit mask),
including the i^omega_L_mult(ms) storage phase. Hermitian for every weight.
"""
function dense_string(ms_int::Integer, nf::Int, gammas::Vector{Matrix{ComplexF64}})
    M = Matrix{ComplexF64}(LinearAlgebra.I, 2^nf, 2^nf)
    for k in 1:2nf
        if (ms_int >> (k - 1)) & 1 == 1
            M = M * gammas[k]
        end
    end
    w = get_weight(ms_int)
    phase = (1im)^Int(mod(div(w * (w - 1), 2), 2))
    return phase * M
end

"""
    dense_sum(msum, gammas)
Dense matrix of an `AbstractMajoranaSum`.
"""
function dense_sum(msum, gammas::Vector{Matrix{ComplexF64}})
    nf = nfermions(msum)
    M = zeros(ComplexF64, 2^nf, 2^nf)
    for (ms, coeff) in zip(terms(msum), coefficients(msum))
        M .+= coeff .* dense_string(ms, nf, gammas)
    end
    return M
end

const _e0 = ComplexF64[1, 0]
const _e1 = ComplexF64[0, 1]

_dense_fock_nf(f::FockState) = f.is_spinful ? 2 * f.n_sites : f.n_sites
_fock_occupied(f::FockState, mode::Int) = (f.occupied_sites >> (2 * mode - 2)) & 1 == 1

"""
    dense_fock(f::FockState)
Dense state vector of a Fock basis state, in the package's phase convention
(descending-ordered creation product, eta = (-1)^(n(n-1)/2)).
"""
function dense_fock(f::FockState)
    nf = _dense_fock_nf(f)
    vecs = [_fock_occupied(f, s) ? _e1 : _e0 for s in 1:nf]
    n = count(s -> _fock_occupied(f, s), 1:nf)
    eta = (-1)^(div(n * (n - 1), 2) % 2)
    return eta * reduce(kron, vecs)
end

"""
    dense_annihilators(nf, gammas)
Dense annihilation operators a_s = (gamma_s + i gamma'_s)/2 for each mode.
For spinful systems, mode 2s-1 is (site s, up) and mode 2s is (site s, down).
"""
dense_annihilators(nf::Int, gammas::Vector{Matrix{ComplexF64}}) =
    [(gammas[2s-1] + im * gammas[2s]) / 2 for s in 1:nf]

"""
    allfockstates(n_sites; is_spinful=false)
All Fock basis states of the system, for exhaustive overlap checks.
"""
function allfockstates(n_sites::Int; is_spinful::Bool=false)
    if is_spinful
        states = FockState[]
        for m in 0:(4^n_sites-1)
            up = [s for s in 1:n_sites if (m >> (2s - 2)) & 1 == 1]
            dn = [s for s in 1:n_sites if (m >> (2s - 1)) & 1 == 1]
            push!(states, FockState(n_sites, up, dn))
        end
        return states
    else
        return [FockState(n_sites, [s for s in 1:n_sites if (m >> (s - 1)) & 1 == 1]) for m in 0:(2^n_sites-1)]
    end
end

"""
    random_string(nf)
Random Majorana string (as an integer, any weight) on nf modes.
"""
function random_string(nf::Int)
    TT = getinttype(nf)
    # mask built bit-by-bit: TT(1) << 2nf overflows to 0 when 2nf == bitwidth(TT)
    mask = zero(TT)
    for k in 1:2nf
        mask |= TT(1) << (k - 1)
    end
    return rand(TT) & mask
end

"""
    random_even_string(nf)
Random non-identity even-weight Majorana string (as an integer) on nf modes.
"""
function random_even_string(nf::Int)
    g = random_string(nf)
    while get_weight(g) % 2 != 0 || g == 0
        g = random_string(nf)
    end
    return g
end

"""
    random_even_msum(nf, nterms; complexcoeffs=false, is_spinful=false)
Random spinless MajoranaSum with `nterms` even-weight strings and randn coefficients.
"""
function random_even_msum(nf::Int, nterms::Int; complexcoeffs::Bool=false, is_spinful::Bool=false)
    CT = complexcoeffs ? ComplexF64 : Float64
    n_sites = is_spinful ? nf ÷ 2 : nf
    msum = MajoranaSum(CT, n_sites, is_spinful)
    for _ in 1:nterms
        c = complexcoeffs ? randn(ComplexF64) : randn()
        PauliPropagation.PropagationBase.add!(msum, random_even_string(nf), c)
    end
    return msum
end
