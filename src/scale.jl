"""
    scaling_factor(A, b; atol = 1e-8)

Return scaling factor used to rescale constraints in matrix format (row-wise).
Arguments are assumed to be of the form `A ∝ b` where `∝` are usually `=` or `≤`
in the context of an optimization problem. To set the cut-off for when a coefficient
is considered to be zero, change `atol`.

See also [`check_scaling`](@ref).
"""
function scaling_factor(A, b; atol = 1e-8)
    mat = [A b]
    max_coeff_range, _ = check_scaling(mat; atol)

    lb = -round(max_coeff_range, RoundUp) / 2.0
    ub = round(max_coeff_range, RoundUp) / 2.0

    sfs = zeros(size(mat, 1))
    for (j, row) in enumerate(eachrow(mat))
        llv, luv = extrema(log10 ∘ abs, row)
        if lb <= llv && luv <= ub
            sf = 0.0
        else # scale to upper bound
            sf = ub - luv
        end
        sfs[j] = 10^sf
    end

    return sfs
end

"""
    _get_exponent_or_cutoff(x; atol = 1e-8)

Return `log₁₀(|x|)` or `x` if `|x| ≤ atol`.
"""
_get_exponent_or_cutoff(x; atol = 1e-8) = abs(x) > atol && log10(abs(x))

"""
    check_scaling(mat; atol = 1e-8)

Return the best case (if rescaled using [`scaling_factor`](@ref)) and current
scaling of a matrix `mat`. Scaling is defined as the largest difference between
the exponent of the smallest and largest value in each row of a matrix. To set
the cut-off for when a coefficient is considered to be zero, change `atol`.
"""
function check_scaling(mat; atol = 1e-8)
    best_case = maximum(
        maximum(_get_exponent_or_cutoff.(mat; atol), dims = 2)[:] -
        minimum(_get_exponent_or_cutoff.(mat; atol), dims = 2)[:],
    )
    rlb, rub = log10.(extrema(filter(x -> x > atol, abs.(mat))))
    worst_case = rub - rlb
    return best_case, worst_case
end
