"""
    $(TYPEDSIGNATURES)

Return the analytic derivative of the optimality conditions with respect to 
the parameters.
"""
function derivative_of_enzyme_equality(gm::GeckoModel, rid_enzyme)
    E_components, kcat_stoich_idx = _build_equality_enzyme_constraints(gm, rid_enzyme)

    row_col_stoich_tidx = zip(
        E_components.row_idxs,
        E_components.col_idxs,
        [stoich for (stoich, _) in kcat_stoich_idx],
        [idx for (_, idx) in kcat_stoich_idx],
    )

    function db1(x, ν, λ, θ, reg)
        deriv_block = zeros(length(x), length(rid_enzyme))
        for j = 1:length(x)
            for (_r, c, s, i) in row_col_stoich_tidx
                r = _r + length(ν) - n_genes(gm)
                if c == j
                    deriv_block[j, i] -= -ν[r] * s / θ[i]^2 # ν is negative by JuMP assumption   
                end
            end
        end
        deriv_block
    end

    function db2(x, ν, λ, θ, reg)
        deriv_block = zeros(n_genes(gm), length(rid_enzyme))
        for j = 1:n_genes(gm)
            for (r, c, s, i) in row_col_stoich_tidx
                if r == j
                    deriv_block[j, i] -= x[c] * s / θ[i]^2
                end
            end
        end
        [
            zeros(length(ν) - n_genes(gm), length(rid_enzyme))
            deriv_block
        ]
    end

    num_ineq_cons = n_genes(gm) * 2 + n_reactions(gm) * 2 + n_coupling_constraints(gm) * 2

    return (x, ν, λ, θ, reg) -> [
        db1(x, ν, λ, θ, reg)
        db2(x, ν, λ, θ, reg)
        zeros(num_ineq_cons, length(rid_enzyme))
    ]
end

