function truncate_hybrid_string(hs::sT, coeff, nq, nf, maxW_qubit, maxW_fermions, maxsingles_fermion, qubits_filter, fermions_filter, fermions_max_single_filter) where {sT<:Integer}
    if maxW_qubit < Inf && PauliPropagation.truncateweight(hs & qubits_filter, maxW_qubit)
        #println("Truncating hybrid string due to qubit weight: $(hs | qubits_filter)")
        return true
    end
    if maxW_fermions < Inf && truncatemajoranaweight(hs & fermions_filter, maxW_fermions)
        return true
    end
    if maxsingles_fermion < Inf && (compute_max_single(hs & fermions_filter, nf, fermions_max_single_filter) > maxsingles_fermion)
        return true
    end



    return false
end