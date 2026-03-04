include("yao_utils.jl")

function circ_to_yao(nq, fg_circ::Vector{FermionicGate}, thetas; mult_pref = 2.)
    @assert length(fg_circ) == length(thetas)
    circ = chain(nq)
    for (gate, theta) in zip(fg_circ, thetas)
        if gate.symbol == :hop 
            pp_rot1 = PauliRotation([:X, :X], gate.sites)
            pp_rot2 = PauliRotation([:Y, :Y], gate.sites)
            elem_to_yao!(circ, pp_rot1, -0.5 * theta * mult_pref)
            elem_to_yao!(circ, pp_rot2, -0.5 * theta * mult_pref)
        elseif gate.symbol == :nn 
            pp_rot = PauliRotation([:Z, :Z], gate.sites)
            pp_rot_single1 = PauliRotation([:Z], gate.sites[1])
            pp_rot_single2 = PauliRotation([:Z], gate.sites[2])
            elem_to_yao!(circ, pp_rot, 0.25 * theta * mult_pref)
            elem_to_yao!(circ, pp_rot_single1, -0.25 * theta * mult_pref)
            elem_to_yao!(circ, pp_rot_single2, -0.25 * theta * mult_pref)
        else 
            error("Unsupported gate symbol: $(gate.symbol)")
        end
    end
    return circ
end