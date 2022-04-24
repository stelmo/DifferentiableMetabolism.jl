"""
    make_analytic_derivatives(diffmodel::DifferentiableModel)

Creates analytic derivative functions of the `diffmodel` using symbolic
variables.
"""
function make_analytic_derivatives(diffmodel::DifferentiableModel)
    dm = diffmodel # convenience

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
    #     # dm.Q(θ) * x + dm.c(θ) - dm.E(θ)' * ν - dm.M(θ)' * λ
    #     dm.E(θ) * x - dm.d(θ)
    #     λ .* (dm.M(θ) * x - dm.h(θ))
    # ]
    # sparse_F(sx, sν, sλ, sθ)

    # rf = sparse_F(sx, sν, sλ, sθ)
    # sj = Symbolics.jacobian(rf, [sx, sν, sλ, sθ])

    sparse_F(z) = [
        dm.Q(z[θidxs]) * z[xidxs] + dm.c(z[θidxs]) - dm.E(z[θidxs])' * z[νidxs] -
        dm.M(z[θidxs])' * z[λidxs]
        dm.E(z[θidxs]) * z[xidxs] - dm.d(z[θidxs])
        spdiagm(z[λidxs]) * (dm.M(z[θidxs]) * z[xidxs] - dm.h(z[θidxs]))
    ]
    sz = [sx; sν; sλ; sθ]
    sparse_F(sz)

    sj = sparse(Symbolics.jacobian(sparse_F(sz), sz)) #TODO only get jacobian of theta

    f_expr = build_function(sj, [sx; sν; sλ; sθ])
    myf = eval(first(f_expr))
    nrows = length(sx) + length(sν) + length(sλ)

    dm.analytic_par_derivs =
        (x, ν, λ, θ) -> reshape(myf([x; ν; λ; θ])[(nrows^2+1):end], nrows, length(sθ))

    dm.analytic_var_derivs =
        (x, ν, λ, θ) -> [
            dm.Q(θ) -dm.E(θ)' -dm.M(θ)'
            dm.E(θ) spzeros(size(dm.E(θ), 1), length(ν)) spzeros(size(dm.E(θ), 1), length(λ))
            spdiagm(λ)*dm.M(θ) spzeros(size(dm.M(θ), 1), length(ν)) spdiagm(dm.M(θ) * x - dm.h(θ))
        ]
    return nothing
end

