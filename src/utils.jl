"""
$(TYPEDSIGNATURES)

Remove linearly dependent rows in `A` deduced by QR decomposition. After QR decomposition, 
all elements are rounded to `digits` after the decimal. Rows of all zeros
(absolute value of the each element in row ≤ `atol`) are removed.
"""
function remove_lin_dep_rows_QR(A; atol = 1e-8)
    R = qr(Array(A)).R

    return R[any(x -> abs(x) > atol, R, dims = 2)[:], :] # remove rows of all zero

end

"""
$(TYPEDSIGNATURES)

Remove linearly dependent rows in `A` deduced by reducing the matrix to row
echelon form. Adjust sensitivity to numerical issues with `ϵ`. After row echelon
reduction, all elements are rounded to `digits` after the decimal. Rows of all zeros
(absolute value of the each element in row ≤ `atol`) are removed. Beware, this
operation is expensive for very large matrices.
"""
function remove_lin_dep_rows(A; ϵ = 1e-8, atol = 1e-8, digits = 16)
    #TODO this method is suboptimal and can be improved for numerical stability
    #TODO improve RowEchelon, SVD does not work due to column reordering
    rA = round.(RowEchelon.rref!(copy(Array(A)), ϵ); digits)

    idxs = Int[]
    for i = 1:size(rA, 1)
        if !all(abs.(rA[i, :]) .<= atol) # remove rows of all zero
            push!(idxs, i)
        end
    end

    return rA[idxs, :]
end

export remove_lin_dep_rows
