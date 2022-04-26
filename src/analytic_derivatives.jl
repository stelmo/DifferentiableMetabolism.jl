"""
    make_symbolic_derivatives(diffmodel::DifferentiableModel)

Creates analytic derivative functions of the `diffmodel` internals, using
symbolic variables. Note, this function can take some time to construct the
derivatives, but substantially speeds up repeated calls to
[`differentiate`](@ref).
"""
function make_symbolic_derivatives(diffmodel::DifferentiableModel)
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

    sparse_F(x, ν, λ, θ) = [
        dm.Q(θ) * x + dm.c(θ) - _transpose(dm.E(θ)) * ν - _transpose(dm.M(θ)) * λ
        dm.E(θ) * x - dm.d(θ)
        λ .* (dm.M(θ) * x - dm.h(θ))
    ] #TODO remove transpose once #575 on Symbolics gets fixed
    # sparse_F(sx, sν, sλ, sθ)

    sj = Symbolics.jacobian(Symbolics.scalarize.(sparse_F(sx, sν, sλ, sθ)), sθ)

    f_expr = build_function(sj, sx, sν, sλ, sθ)
    myf = eval(first(f_expr))

    dm.analytic_par_derivs = (x, ν, λ, θ) -> myf(x, ν, λ, θ)

    dm.analytic_var_derivs =
        (x, ν, λ, θ) -> [
            dm.Q(θ) -dm.E(θ)' -dm.M(θ)'
            dm.E(θ) spzeros(size(dm.E(θ), 1), length(ν)) spzeros(size(dm.E(θ), 1), length(λ))
            spdiagm(λ)*dm.M(θ) spzeros(size(dm.M(θ), 1), length(ν)) spdiagm(dm.M(θ) * x - dm.h(θ))
        ]
    return nothing
end

function _transpose(mat)
    I, J, V = findnz(mat)
    sparse(J, I, V)
end
