"""
    make_symbolic_derivatives(diffmodel::DifferentiableModel)

Creates analytic derivative functions of the `diffmodel` internals, using
symbolic variables. Note, this function can take some time to construct the
derivatives, but substantially speeds up repeated calls to
[`differentiate`](@ref).
"""
function make_symbolic_derivatives(diffmodel::DifferentiableModel)
    diffmodel.analytic_par_derivs = make_symbolic_param_derivative(diffmodel)
    diffmodel.analytic_var_derivs = make_symbolic_var_derivative(diffmodel)
    return nothing
end

function _transpose(mat)
    I, J, V = findnz(mat)
    sparse(J, I, V)
end

function make_symbolic_param_derivative(dm::DifferentiableModel)
    xidxs = 1:length(dm.c(dm.θ))
    νidxs = last(xidxs) .+ (1:length(dm.d(dm.θ)))
    λidxs = last(νidxs) .+ (1:length(dm.h(dm.θ)))
    θidxs = last(λidxs) .+ (1:length(dm.θ))

    Symbolics.@variables begin
        sx[1:length(xidxs)]
        sν[1:length(νidxs)]
        sλ[1:length(λidxs)]
        sθ[1:length(θidxs)]
    end

    # sparse_F(x, ν, λ, θ) = [
    #     dm.Q(θ) * x + dm.c(θ) - _transpose(dm.E(θ)) * ν - _transpose(dm.M(θ)) * λ
    #     dm.E(θ) * x - dm.d(θ)
    #     λ .* (dm.M(θ) * x - dm.h(θ))
    # ] #TODO remove transpose once #575 on Symbolics gets fixed
    # # sparse_F(sx, sν, sλ, sθ)

    # sj = Symbolics.sparsejacobian(Symbolics.scalarize.(sparse_F(sx, sν, sλ, sθ)), sθ)

    # f_expr = build_function(sj, sx, sν, sλ, sθ)
    # myf = eval(first(f_expr))

    # (x, ν, λ, θ) -> myf(x, ν, λ, θ)

    sparse_F(z) = [
        dm.Q(z[θidxs]) * z[xidxs] + dm.c(z[θidxs]) - dm.E(z[θidxs])' * z[νidxs] - dm.M(z[θidxs])' * z[λidxs]
        dm.E(z[θidxs]) * z[xidxs] - dm.d(z[θidxs])
        z[λidxs] .* (dm.M(z[θidxs]) * z[xidxs] - dm.h(z[θidxs]))
    ]
    sz = [sx; sν; sλ; sθ]
    # sparse_F(sz)

    #TODO only get jacobian of theta
    sj = sparse(Symbolics.jacobian(sparse_F(sz), sz)[:, end-(length(sθ)-1):end])
    (nr, nc) = size(sj)

    f_expr = build_function(sj, [sx; sν; sλ; sθ])
    myf = eval(first(f_expr))

    (x, ν, λ, θ) -> reshape(myf([x; ν; λ; θ]), nr, nc)
end

make_symbolic_var_derivative(dm::DifferentiableModel) = 
    (x, ν, λ, θ) -> [
        dm.Q(θ) -dm.E(θ)' -dm.M(θ)'
        dm.E(θ) spzeros(size(dm.E(θ), 1), length(ν)) spzeros(size(dm.E(θ), 1), length(λ))
        spdiagm(λ)*dm.M(θ) spzeros(size(dm.M(θ), 1), length(ν)) spdiagm(dm.M(θ) * x - dm.h(θ))
    ]

# function _tranpose_mul_vec(A, b)
#     rows = rowvals(A)
#     vals = nonzeros(A)
#     _, n = size(A)
#     for j = 1:n
#        for i in nzrange(A, j)
#           row = rows[i]
#           val = vals[i]
#           # perform sparse wizardry...
#        end
#     end
# end