module MajoranaPropagation

using PauliPropagation
using PauliPropagation.PropagationBase
import PauliPropagation.PropagationBase: propagate, propagate!

include("MajoranaDataTypes.jl")
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
    add!
include("MajoranaAlgebra.jl")
export
    fock_mask,
    overlap_with_fock,
    overlap_with_fock_spinful,
    ms_mult,
    commutator,
    commutes,
    norm,
    omega_mult,
    omega_L_mult

include("propagationcache.jl")
export MajoranaPropagationCache

include("gates.jl")
export
    MajoranaRotation,
    FermionicGate,
    getnewmajoranastring,
    MajoranaRotation,
    countparameters,
    propagate,
    propagate!

include("truncations.jl")
export
    create_unpaired_mask,
    create_doublons_filters,
    compute_unpaired,
    compute_doublons,
    truncatemajoranaweight

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
end