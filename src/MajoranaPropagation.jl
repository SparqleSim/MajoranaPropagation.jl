module MajoranaPropagation

using PauliPropagation
#using ThreadPools
using Base.Threads
import PauliPropagation: set!, coefftype, countparameters, similar, propagate!, applytoall!, applymergetruncate!, checktruncationonall!, mergeandempty!, wrapcoefficients, truncatemincoeff, truncatefrequency, truncatesins, _checkfreqandsinfields, _checkcircandthetas, _promotecircandthetas

include("MajoranaAlgebra.jl")
export
    MajoranaSum,
    MajoranaString,
    nfermions,
    set!,
    length,
    get_weight,
    coefftype,
    similar,
    iterate,
    fock_filter,
    overlap_with_fock,
    overlap_with_fock_spinful,
    getinttype,
    ms_mult,
    add!,
    commutator,
    commutes,
    norm,
    omega_mult,
    omega_L_mult

include("truncations.jl")
export
    checktruncationonall!,
    create_max_single_filter,
    create_doublons_filters,
    compute_max_single,
    compute_doublons,
    truncatemajoranaweight,
    checktruncationonall!

include("gates.jl")
export
    MajoranaRotation,
    FermionicGate,
    applytoall!,
    getnewmajoranastring,
    MajoranaRotation,
    countparameters,
    applymergetruncate!,
    mergeandempty!,
    empty!


include("circuits.jl")
export
    hubbard_circ_fermionic_sites,
    hubbard_circ_fermionic_sites_single_layer,
    fermionic_hubbard_circ_fermionic_sites_single_layer,
    hubbard_circ_fermionic_sites_second_order,
    fermionic_hubbard_circ_fermionic_sites_second_order_single_layer

include("MajoranaFrequencyTracker.jl")
export
    MajoranaFrequencyTracker,
    wrapcoefficients,
    reset_tracker!

include("Constructors.jl")

include("multidict.jl")
export
    MajoranaSumMulti, 
    propagate!,
    show_stats,
    applymergetruncate!,
    mergeandempty!,
    applytoall!
end