const Chol = if isdefined(Base, :SparseArrays)
    Base.SparseArrays.CHOLMOD
else
    Base.SparseMatrix.CHOLMOD
end

colptr0(A::Chol.Sparse, col) = unsafe_load(unsafe_load(A.p).p, col)
rowval(A::Chol.Sparse, j) = unsafe_load(unsafe_load(A.p).i, j) + 1
nzval(A::Chol.Sparse, j) = unsafe_load(unsafe_load(A.p).x, j)

function calcx(E, m, n, k)
    CSE = Chol.Sparse(E)
    tmp = CSE'*CSE

    # subtract spdiagm( diag(tmp) ) from tmp in-place by setting diagonals to 0
    tmp_ptr = unsafe_load(tmp.p)
    for col in 1:n
        for j in colptr0(tmp, col) + 1 : colptr0(tmp, col+1)
            if rowval(tmp, j) == col
                unsafe_store!(tmp_ptr.x, 0, j)
            end
        end
    end

    R = CSE * tmp
    s = zeros(Int, m)

    @inbounds for col in 1:n
        for j in colptr0(R, col) + 1 : colptr0(R, col+1)
            if nzval(R, j) == 2
                s[rowval(R, j)] += 1
            end
        end
    end
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
    #E_colptr = E.colptr
    #E_rowval = E.rowval
    #E_nzval  = E.nzval

    #
    tic()
    m,n = size(E)
    x, xc = calcx(E, m, n, k)
    while sum(xc) != sum( any(E,2) )
        # set elements of E in rows where x is true to 0, E[find(x), :] = 0
        xcrows = find(xc)
        Ekeep = E[xcrows, :]
        E = SparseMatrixCSC(m, n, Ekeep.colptr, xcrows[Ekeep.rowval], Ekeep.nzval)
        #num_deleted = 0
        #@inbounds for col in 1:size(E, 2)
        #    for k in E_colptr[col] + num_deleted : E_colptr[col+1]-1
        #        if x[E_rowval[k]]
        #            num_deleted += 1
        #        elseif num_deleted > 0
        #            E_rowval[k - num_deleted] = E_rowval[k]
        #            E_nzval[k - num_deleted] = E_nzval[k]
        #        end
        #    end
        #    E_colptr[col+1] -= num_deleted
        #end
        #resize!(E_rowval, length(E_rowval) - num_deleted)
        #resize!(E_nzval, length(E_nzval) - num_deleted)
        #E[find(x), :] = 0
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
