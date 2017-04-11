function simpler_readdlm(f, delim, eltyp)
    results = []
    linelength::Int = 0
    open(f) do handle
        for line in eachline(handle)
            curline = Vector{eltyp}(0)
            for elem in split(line, delim)
                push!(curline, parse(eltyp, elem))
            end
            if linelength == 0
                linelength = length(curline)
            elseif linelength != length(curline)
                error("unexpected line length $linelength")
            end
            push!(results, curline)
        end
    end
    output = Matrix{eltyp}(length(results), linelength)
    @inbounds for j = 1:linelength, i = 1:length(results)
        output[i,j] = results[i][j]
    end
    return output
end

function calcx(E, m, n, k)
    tmp = E.'*E

    # subtract spdiagm( diag(tmp) ) from tmp in-place by setting diagonals to 0
    # hoist field access (shouldn't be necessary on julia >= 0.5)
    tmp_colptr = tmp.colptr
    tmp_rowval = tmp.rowval
    tmp_nzval  = tmp.nzval
    @inbounds for col in 1:size(tmp, 2)
        k1 = tmp_colptr[col]
        k2 = tmp_colptr[col+1]-1
        (k1 > k2) && continue # empty column
        k1 = searchsortedfirst(tmp_rowval, col, k1, k2, Base.Order.Forward)
        if k1 <= k2 && tmp_rowval[k1] == col
            tmp_nzval[k1] = 0
        end
        #for k in tmp_colptr[col] : tmp_colptr[col+1]-1
        #    if tmp_rowval[k] == col
        #        tmp_nzval[k] = 0
        #    elseif tmp_rowval[k] > col
        #        break
        #    end
        #end
    end

    R = E * tmp
    # set elements where E[i,j]==2 to 1, and otherwise to 0 in-place
    # hoist field access (shouldn't be necessary on julia >= 0.5)
    #R_colptr = R.colptr
    #R_rowval = R.rowval
    #R_nzval  = R.nzval
    #@inbounds for col in 1:size(R, 2)
    #    for k in R_colptr[col] : R_colptr[col+1]-1
    #        if R_nzval[k] == 2
    #            R_nzval[k] = 1
    #        else
    #            R_nzval[k] = 0
    #        end
    #    end
    #end
    #s = sum(R, 2)

    #R = E * ( tmp - spdiagm( diag(tmp) ) )
    r,c,v = findnz(R)
    id = v.==2
    A = sparse( r[id], c[id], 1, m, n)
    s = sum(A, 2)
    x = s .< (k-2)

    return (x, !x)
end

function ktruss(inc_mtx_file, k)
    if !isfile( inc_mtx_file )
        println("unable to open input file")
        return (-1)
    end

    # load input data       
    t_read_inc=@elapsed ii = readdlm( inc_mtx_file, '\t', Int64)
    println("incidence matrix read time : ", t_read_inc)

    t_create_inc=@elapsed E = sparse( ii[:,1], ii[:,2], ii[:,3] )
    println("sparse adj. matrix creation time : ", t_create_inc)

    # hoist field access (shouldn't be necessary on julia >= 0.5)
    E_colptr = E.colptr
    E_rowval = E.rowval
    E_nzval  = E.nzval

    #
    tic()
    m,n = size(E)
    x, xc = calcx(E, m, n, k)
    while sum(xc) != sum( any(E,2) )
        # set elements of E in rows where x is true to 0, E[find(x), :] = 0
        #xrows = find(x)
        #@inbounds for col in 1:size(E, 2)
        #    for k in E_colptr[col] : E_colptr[col+1]-1
        #        if E_rowval[k] in xrows
        #            E_nzval[k] = 0
        #        end
        #    end
        #end
        E[find(x), :] = 0
        x, xc = calcx(E, m, n, k)
    end
    
    return E
end



#######################################################
# Graph Challenge Benchmark
# Architect : Dr. Jeremy Kepner (kepner@ll.mit.edu)
# Developer : Dr. Siddharth Samsi (ssamsi@mit.edu)
#
# MIT
########################################################
# (c) <2017> Massachusetts Institute of Technology
########################################################
