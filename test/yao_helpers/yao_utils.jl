using Yao 

function circ_to_yao(nq, pp_circ, thetas)
    @assert length(pp_circ) == length(thetas)
    circ = chain(nq)
    for (gate, theta) in zip(pp_circ, thetas)
        elem_to_yao!(circ, gate, theta)
    end
    return circ
end

function elem_to_yao!(circuit, elem, theta)
    if elem.symbols == [:X, :X]
        push!(circuit, cnot(elem.qinds[2], elem.qinds[1]))
        push!(circuit, Yao.put(elem.qinds[2]=>Yao.Rx(theta)))
        push!(circuit, cnot(elem.qinds[2], elem.qinds[1]))
    elseif elem.symbols == [:Y, :Y]
        push!(circuit, Yao.put(elem.qinds[1]=>Yao.Rz(-pi/2)))
        push!(circuit, Yao.put(elem.qinds[2]=>Yao.Rz(-pi/2)))
        push!(circuit, cnot(elem.qinds[2], elem.qinds[1]))
        push!(circuit, Yao.put(elem.qinds[2]=>Yao.Rx(theta)))
        push!(circuit, cnot(elem.qinds[2], elem.qinds[1]))
        push!(circuit, Yao.put(elem.qinds[1]=>Yao.Rz(pi/2)))
        push!(circuit, Yao.put(elem.qinds[2]=>Yao.Rz(pi/2)))
    elseif elem.symbols == [:Z, :Z]
        push!(circuit, cnot(elem.qinds[1], elem.qinds[2]))
        push!(circuit, Yao.put(elem.qinds[2]=>Yao.Rz(theta)))
        push!(circuit, cnot(elem.qinds[1], elem.qinds[2]))
    elseif elem.symbols == [:X]
        push!(circuit, Yao.put(elem.qinds[1]=>Yao.Rx(theta)))
    elseif elem.symbols == [:Y]
        push!(circuit, Yao.put(elem.qinds[1]=>Yao.Ry(theta)))
    elseif elem.symbols == [:Z]
        push!(circuit, Yao.put(elem.qinds[1]=>Yao.Rz(theta)))
    else
        error("Unsupported gate symbols: $(elem.symbols)")
    end 
    return 
end
