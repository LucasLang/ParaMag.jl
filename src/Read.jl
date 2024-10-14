"""
Read parameters (in atomic units).
method: e.g. "CASSCF" or "NEVPT2"
TO DO: extend for f elements and for SOC parameter.
"""
function read_AILFT_params_ORCA(outfile::String, method::String)
    nel = parse_int(outfile, ["nel"], 0, 3)
    norb = parse_int(outfile, ["norb"], 0, 3)
    l = (norb-1)÷2
    if norb == 5
        hLFT = Matrix{Float64}(undef, norb, norb)
        for row in 1:norb
            for col in 1:norb
                hLFT[row,col] = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Orbital"], row, col)
            end
        end
        perm = [4,2,1,3,5]    # change order from 0,1,-1,2,-2 to 2,1,0,-1,-2 (=x2-y2,xz,z2,yz,xy) 
        hLFT = hLFT[perm, perm]
        F0 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 2, 4)
        F2 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 3, 2)
        F4 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 4, 2)
        F = Dict(0 => F0, 2 => F2/49, 4 => F4/441)
        zeta = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "ZETA_D"], 0, 2)/219474.63  # convert from cm-1 to Hartree
    end
    if norb == 7
        hLFT = Matrix{Float64}(undef, norb, norb)
        for row in 1:norb
            for col in 1:norb
                hLFT[row,col] = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Orbital"], row, col)
            end
        end
        perm = [6,4,2,1,3,5,7]    # change order from 0,1,-1,2,-2,3,-3 to 3,2,1,0,-1,-2,-3
        hLFT = hLFT[perm, perm]
        F0 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 2, 2)
        F2 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 3, 2)
        F4 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 4, 2)
        F6 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 5, 2)
        F = Dict(0 => F0, 2 => F2/15/15, 4 => F4/33/33, 6 => (5/429)^2 * F6)
        zeta = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "ZETA_F"], 0, 2)/219474.63  # convert from cm-1 to Hartree
    end
    return LFTParam(nel, norb, l, hLFT, F, zeta)
end

function read_AILFT_params_ORCA6(outfile::String, method::String)
    nel = parse_int(outfile, ["nel"], 0, 3)
    norb = parse_int(outfile, ["norb"], 0, 3)
    l = (norb-1)÷2
    if norb == 5
        hLFT = Matrix{Float64}(undef, norb, norb)
        for row in 1:norb
            for col in 1:norb
                hLFT[row,col] = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Orbital"], row, col)
            end
        end
        perm = [4,2,1,3,5]    # change order from 0,1,-1,2,-2 to 2,1,0,-1,-2 (=x2-y2,xz,z2,yz,xy) 
        hLFT = hLFT[perm, perm]
<<<<<<< HEAD
        F2 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 2, 2)
        F4 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 3, 2)
        F = Dict(0 => 0, 2 => F2/49, 4 => F4/441)
=======
        if method == "CASSCF"
            F0 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 2, 4)
            F2 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 3, 2)
            F4 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 4, 2)
        else
            F0 = 0.0
            F2 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 2, 2)
            F4 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 3, 2)
        end
        F = Dict(0 => F0, 2 => F2/49, 4 => F4/441)
        println("F ", F0," ", F2," ", F4)
>>>>>>> calcshift
        zeta = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "ZETA_D"], 0, 2)/219474.63  # convert from cm-1 to Hartree
    end
    if norb == 7
        hLFT = Matrix{Float64}(undef, norb, norb)
        for row in 1:norb
            for col in 1:norb
                hLFT[row,col] = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Orbital"], row, col)
            end
        end
        
        perm = [6,4,2,1,3,5,7]   # change order from 0,1,-1,2,-2,3,-3 to 3,2,1,0,-1,-2,-3
    
        hLFT = hLFT[perm, perm]
<<<<<<< HEAD
=======

        if method != "CASSCF"
            println("WARNING !!!!!!")
        end
        
>>>>>>> calcshift
        F0 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 2, 2)
        F2 = 0
        try
            F2 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 3, 2)
        catch
        end
        F4 = 0
        F6 = 0
        try
            F4 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 4, 2)
            F6 = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "Slater-Condon"], 5, 2)
        catch
        end
        F = Dict(0 => F0, 2 => F2/15/15, 4 => F4/33/33, 6 => (5/429)^2 * F6)
        zeta = parse_float(outfile, ["AILFT MATRIX ELEMENTS ($method)", "ZETA_F"], 0, 2)/219474.63  # convert from cm-1 to Hartree
    end
    return LFTParam(nel, norb, l, hLFT, F, zeta)
end

function is_float(element)
    try
        parse(Float64, element)
        return true
    catch
        return false
    end
end


function read_Aiso(filecalcnmr::String)

    file = readlines(filecalcnmr)
    Aiso = Vector{Float64}(undef, 0)
    for (i, line) in enumerate(file)

        if occursin("A(iso)", line) 
            splitline = split(line, "=")
            push!(Aiso, parse(Float64, splitline[end]))
        end
    end

    Aiso = reshape(Aiso ,(1,length(Aiso)))
    return Aiso
end


function read_HFCmatrix(filecalcnmr::String, natoms::Int)
    #works for orca6

    file = readlines(filecalcnmr)

    HFC = []
    AFC = []
    count = 0
    for (i, line) in enumerate(file)

        if count > natoms
            count = 0
            HFC = []
            AFC = []
        end

        if occursin("Raw HFC matrix (all values in MHz)", line) && count < natoms
            count += 1
            v_single = zeros(Float64, 3, 3)
            for (ij,j) in enumerate(2:1:4)
                splitline = split(file[i + j], "         ")
                lista = [parse(Float64, splitline[l]) for l in 2:length(splitline) if is_float(splitline[l])]
                for jj in eachindex(lista)
                    v_single[ij,jj] = lista[jj]
                end
            end
            push!(HFC, v_single)
        end

        if occursin("A(FC)", line) && count < natoms+1
            splitline = split(line, "    ")
            push!(AFC, [parse(Float64, splitline[j]) for j in 2:length(splitline) if is_float(splitline[j])])
        end

    end

    Aorb = []
    for i in eachindex(AFC)
        push!(Aorb, HFC[i] .- Diagonal(AFC[i]))
    end 

    return Aorb
end


<<<<<<< HEAD
function read_effectiveH(filecalcnmr::String, theory::String, Dflag::Bool=false, gflag::Bool=false)
=======
function read_effectiveH(filecalcnmr::String, theory::String, Smult::Int, Dflag::Bool=false, gflag::Bool=false)
>>>>>>> calcshift
    #works for orca6

    file = readlines(filecalcnmr)

    D_matrices = []
    g_matrices = []

    for (i, line) in enumerate(file)

        if Dflag
<<<<<<< HEAD
            if occursin("Raw matrix (cm-1)", line) && occursin("ZERO-FIELD SPLITTING", file[i-19]) && occursin("EFFECTIVE HAMILTONIAN", file[i-18])  #depends on S
=======
            if occursin("Raw matrix (cm-1)", line) && occursin("ZERO-FIELD SPLITTING", file[i-((2*Smult)+11)]) && occursin("EFFECTIVE HAMILTONIAN", file[i-((2*Smult)+10)])  #depends on S
>>>>>>> calcshift
                matrix = zeros(Float64, 3, 3)
                for (ij,j) in enumerate(1:1:3)
                    splitline = split(file[i + j], "  ")
                    lista = [parse(Float64, splitline[l]) for l in 2:length(splitline) if is_float(splitline[l])]
                    for jj in eachindex(lista)
                        matrix[ij,jj] = lista[jj]
                    end
                end
                push!(D_matrices, matrix)
            end
        end

        if gflag
            if occursin("g-matrix:", line) && occursin("ELECTRONIC G-MATRIX FROM EFFECTIVE HAMILTONIAN", file[i-4])
                matrix = zeros(Float64, 3, 3)
                for (ij,j) in enumerate(1:1:3)
                    splitline = split(file[i + j], "  ")
                    lista = [parse(Float64, splitline[l]) for l in 2:length(splitline) if is_float(splitline[l])]
                    for jj in eachindex(lista)
                        matrix[ij,jj] = lista[jj]
                    end
                end
                push!(g_matrices, matrix)
            end
        end
    end

    returned = []
    if theory == "NEVPT2"
        if Dflag
            push!(returned, D_matrices[2])
        end
        if gflag
            push!(returned, g_matrices[2])
        end
    elseif theory == "CASSCF"
        if Dflag
            push!(returned, D_matrices[1])
        end
        if gflag
            push!(returned, g_matrices[1])
        end
    else
        if Dflag
            push!(returned, D_matrices)
        end
        if gflag
            push!(returned, g_matrices)
        end
    end
    return returned
end


function read_DFT(filecalcnmr::String, Dflag::Bool=false, gflag::Bool=false)
    #works for orca6

    file = readlines(filecalcnmr)

    D_matrix = Matrix{Float64}(undef, 0, 0)
    g_matrix = Matrix{Float64}(undef, 0, 0)

    for (i, line) in enumerate(file)

        if Dflag
            if occursin("raw-matrix :", line) && occursin("ZERO-FIELD-SPLITTING TENSOR", file[i-3])
                D_matrix = zeros(Float64, 3, 3)
                for (ij,j) in enumerate(1:1:3)
                    splitline = split(file[i + j], "  ")
                    lista = [parse(Float64, splitline[l]) for l in 2:length(splitline) if is_float(splitline[l])]
                    for jj in eachindex(lista)
                        D_matrix[ij,jj] = lista[jj]
                    end
                end
            end
        end

        if gflag
            if occursin("g-matrix:", line) && occursin("ELECTRONIC G-MATRIX", file[i-3])
                g_matrix = zeros(Float64, 3, 3)
                for (ij,j) in enumerate(1:1:3)
                    splitline = split(file[i + j], "  ")
                    lista = [parse(Float64, splitline[l]) for l in 2:length(splitline) if is_float(splitline[l])]
                    for jj in eachindex(lista)
                        g_matrix[ij,jj] = lista[jj]
                    end
                end
            end
        end
    end

    if gflag==false && Dflag
        return D_matrix
    elseif Dflag==false && gflag
        return g_matrix
    elseif gflag && Dflag
        return D_matrix, g_matrix
    end

end


function read_integrals_oneelint(fileint::String, dim::Int)
    file = readlines(fileint)
    integrals = Float64[]
    checkpoint=false
    for (i, line) in enumerate(file)

        if checkpoint
            splitline = split(line, "    ")
            push!(integrals, parse(Float64, splitline[end]))
        end

        if occursin("one-electron integrals", line)
            checkpoint=true
        end

    end
    integrals = reshape(integrals, (dim,dim)) 
    
    perm = [6,4,2,1,3,5,7]    # change order from 0,1,-1,2,-2,3,-3 to 3,2,1,0,-1,-2,-3
    integrals = integrals[perm, perm]

    return integrals
end


function read_integrals_ee(fileint::String, dim::Int)

    file = readlines(fileint)
    integrals = Float64[]
    DIM = dim*dim*dim*dim
    count = 0
    for (i,line) in enumerate(file)
        splitline = split(line, "    ")
        if count < DIM
            push!(integrals, parse(Float64, splitline[end]))
        else
            break
        end
        count+=1
    end
    integrals = reshape(integrals, (dim,dim,dim,dim))
    perm = [6,4,2,1,3,5,7]    # change order from 0,1,-1,2,-2,3,-3 to 3,2,1,0,-1,-2,-3

    # for i in 1:size(integrals, 1)
    #     for j in 1:size(integrals, 2)
    #         matrice = integrals[i,j,:,:]
    #         matrice = matrice[perm, perm]
    #         integrals[i,j,:,:] = matrice
    #     end
    # end

    # for i in 1:size(integrals, 1)
    #     matrice = integrals[i,:,:,:]
    #     matrice = matrice[perm, :, :]
    #     integrals[i,:,:,:] = matrice
    # end

    # integrals = integrals[perm, :, :, :]

    #equivalent to the implementation above
    integrals = integrals[perm, perm, perm, perm]

    return integrals

end


function read_integrals_so_f(fileint::String)   #could be adjuted to be used for d

    file = readlines(fileint)
    Integrals = NTuple{3, Matrix{Float64}}[]
    for (i, line) in enumerate(file)
        integrals = Float64[]
        if occursin("AI-SOC-", line)
            column = Float64[]
            for j in range(2, 16)  #change here for d configurations
                splitline = split(file[i+j], " ")
                row = Float64[]
                for k in eachindex(splitline)
                    if is_float(splitline[k])
                        push!(row, parse(Float64, splitline[k]))
                    end
                end
                if length(row) > 2
                    if integrals == []
                        integrals = row[2:end]
                    else
                        integrals = hcat(integrals, row[2:end])
                    end
                elseif length(row) == 2
                    push!(column, row[end])
                end
            end
            integrals = hcat(integrals', column)
            if occursin("AI-SOC-X integrals (cm-1)", line)
                Integrals = (integrals, Integrals...)
            elseif occursin("AI-SOC-Y integrals (cm-1)", line)
                Integrals = (Integrals..., integrals)
            elseif occursin("AI-SOC-Z integrals (cm-1)", line)
                Integrals = (Integrals..., integrals)
            end
        end
    end

    Integrals = [complex(matrix) for matrix in Integrals]

    perm = [6,4,2,1,3,5,7]    # change order from 0,1,-1,2,-2,3,-3 to 3,2,1,0,-1,-2,-3
    for i in eachindex(Integrals)
        Integrals[i] = Integrals[i][perm, perm]*im
    end

    # for i in eachindex(Integrals)
    #     Integrals[i] = Integrals[i]*im
    # end
    return Integrals

<<<<<<< HEAD
end
=======
end
>>>>>>> calcshift
