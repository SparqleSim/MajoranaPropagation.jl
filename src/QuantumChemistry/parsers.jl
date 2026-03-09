function parse_fcidump(filename; min_coeff=-1.0, remove_identity=true)
    nsites = 0
    for line in eachline(filename)
        if occursin("NORB", line)
            line = filter(!isspace, line)
            data = split(line, "NORB=")
            nsites = parse(Int, split(data[2], ",")[1])
        end
    end

    is_spinful = true
    H = MajoranaSum(Float64, nsites, is_spinful)


    # read file line by line
    found_end = false
    for line in eachline(filename)
        if occursin("&END", line)
            found_end = true
            continue
        end
        if found_end
            data = split(line)
            indices = parse.(Int, data[2:5])
            val = parse(Float64, data[1])

            val = Float64(val)

            #check if 2 body term 
            if 0 in indices
                indices = indices[1:2]
                term = two_body_term(nsites, indices, 1.0)
                H = H + val * term
                #check if 4 body term
            else
                p, q, r, s = indices
                val *= 0.5
                term = four_body_term(nsites, p, q, r, s)
                H = H + val * term
            end
        end
    end

    for (ms, coeff) in H
        if (abs(coeff) < min_coeff) || (ms == 0 && remove_identity)
            delete!(H, ms)
        end
    end
    return H


end
