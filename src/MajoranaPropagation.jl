module MajoranaPropagation

using PauliPropagation
export bricklayertopology, staircasetopology, rectangletopology, getinttype, terms, coefficients
using PauliPropagation.PropagationBase
import PauliPropagation.PropagationBase: propagate, propagate!

include("MajoranaDataTypes.jl")
export
    MajoranaSum,
    MajoranaString,
    nfermions,
    set!,
    add!,
    length,
    get_weight,
    coefftype,
    similar,
    VectorMajoranaSum,
    storage,
    resize!

include("MajoranaAlgebra.jl")
export
    ms_mult,
    commutator,
    commutes,
    norm,
    omega_mult,
    omega_L_mult,
    scalarproduct

include("initial_states.jl")
export
fock_mask,
    FockState,
    overlapwithfock

include("propagationcache.jl")
export
    MajoranaPropagationCache,
    nfermions

include("gates.jl")
export
    MajoranaRotation,
    FermionicRotation,
    getnewmajoranastring,
    MajoranaRotation,
    countparameters,
    propagate,
    propagate!

include("imaginary_gates.jl")
export
    ImaginaryMajoranaRotation,
    ImaginaryFermionicGate

include("propagation.jl")
export
    propagate,
    propagate!

include("truncations.jl")
export
    create_unpaired_mask,
    create_doublons_filters,
    compute_unpaired,
    compute_doublons,
    truncatemajoranaweight,
    truncate!

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

include("QuantumChemistry/QuantumChemistry.jl")
using .QuantumChemistry
end