
function hybridSum(n_qubits::Integer, q_symbols::Vector{Symbol}, q_indices::Vector{Int}, n_fermionic_sites::Integer, f_symb::Symbol, f_sites)
    fermionic_part = MajoranaSum(n_fermionic_sites, Val(f_symb), f_sites)
    TT = getinttype(n_qubits + nfermions(fermionic_part))

    pstr = PauliString(n_qubits + nfermions(fermionic_part), q_symbols, q_indices .+ nfermions(fermionic_part)).term

    hsum_dict = Dict{TT, coefftype(fermionic_part)}()
    for (ms, coeff) in fermionic_part 
        hs = pstr | convert(TT, ms)
        hsum_dict[hs] = coeff
    end
    return hybridSum(n_qubits, n_fermionic_sites, fermionic_part.is_spinful, hsum_dict)
end

function hybridSum(n_qbits::Integer, q_symbols::Vector{Symbol}, q_indices::Vector{Int}, n_fermionic_sites::Integer, is_spinful::Bool)
    n_fermions = is_spinful ? 2 * n_fermionic_sites : n_fermionic_sites
    TT = getinttype(n_qbits + n_fermions)
    pstr = PauliString(n_qbits + n_fermions, q_symbols, q_indices .+ n_fermions).term
    hsum_dict = Dict{TT, Float64}()
    hsum_dict[pstr] = 1.0
    return hybridSum(n_qbits, n_fermionic_sites, is_spinful, hsum_dict)
end