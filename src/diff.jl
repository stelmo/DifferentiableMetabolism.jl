"""
Implicitly differentiate the system and scale.

"""
function differentiate_LP(
    opt,
    θ,
    optimizer;
    use_analytic = false,
    scale_input = true,
    scale_output = true,
    modifications = [],
    dFdθ = nothing,
)
    #: forward pass, solve the optimization problem
    Ef = opt.Ef
    hf = opt.hf
    M = opt.M
    d = opt.d

    E = Ef(θ)
    if scale_input
        !use_analytic &&
            throw(ArgumentError("Can only scale if user supplies parameter derivative.")) # can only scale if using analytic deriv
        E, d = rescale(E, d)
    end

    h = hf(θ)

    opt_model = Model(optimizer)
    set_silent(opt_model)
    x = @variable(opt_model, x[1:size(E, 2)])
    @objective(opt_model, Min, opt.c' * x)
    @constraint(opt_model, eq, E * x .== d)
    @constraint(opt_model, ineq, M * x .<= h)

    # apply the modifications, if any (from COBREXA, only optimizer based mods though)
    for mod in modifications
        mod(nothing, opt_model)
    end

    optimize!(opt_model)

    #: differentiate the optimal solution
    z = value.(opt_model[:x])
    n_vars = length(z)
    ν = dual.(opt_model[:eq])
    n_ν = length(ν)
    λ = dual.(opt_model[:ineq])
    n_λ = length(λ)

    vars = [z; ν; λ]

    if use_analytic
        #: faster
        n_mets = size(E, 1)
        n_rxns = size(E, 2)
        n_cons = size(M, 1)

        A = [
            zeros(n_rxns, n_rxns) -E' -M'
            E zeros(n_mets, n_ν) zeros(n_mets, n_λ)
            diagm(λ)*M zeros(n_cons, n_ν) diagm(M * z - h)
        ]
        isnothing(dFdθ) && throw(
            ArgumentError("Need to supply analytic derivative of parameter functions."),
        )
        B = dFdθ(z, ν, λ, θ)
    else
        #: slower (need to cast to arrays for ForwardDiff)
        F(x, θ) = [
            Array(opt.c) - Ef(θ)' * x[n_vars.+(1:n_ν)] - M' * x[(n_vars+n_ν).+(1:n_λ)]
            Ef(θ) * x[1:n_vars] - Array(d)
            diagm(x[(n_vars+n_ν).+(1:n_λ)]) * (Array(M) * x[1:n_vars] - hf(θ))
        ]
        A = ForwardDiff.jacobian(x -> F(x, θ), vars)
        B = ForwardDiff.jacobian(θ -> F(vars, θ), θ)
    end

    dx = -sparse(A) \ B #! will fail if det(A) = 0
    dx = dx[1:n_vars, :] # only return derivatives of variables, not the duals

    # Scale dx/dy => dlog(x)/dlog(y)
    if scale_output
        scaled_dx = similar(dx)
        for i = 1:size(dx, 1)
            for j = 1:size(dx, 2)
                scaled_dx[i, j] = round(θ[j] / vars[i] * dx[i, j]; digits = 8)
            end
        end
        return value.(x), scaled_dx, objective_value(opt_model)
    else
        return value.(x), dx, objective_value(opt_model)
    end
end

"""
Implicitly differentiate the system and scale.
"""
function differentiate_QP(
    Q,
    c,
    n,
    opt,
    θ,
    optimizer;
    use_analytic = true,
    scale_input = false,
    scale_output = true,
    modifications = [],
    dFdθ = nothing,
)
    Ef = opt.Ef
    hf = opt.hf
    M = opt.M
    d = opt.d

    E = Ef(θ)
    if scale_input
        E, d = rescale(E, d)
    end
    h = hf(θ)

    #: forward pass
    opt_model = Model(optimizer)
    x = @variable(opt_model, x[1:size(E, 2)])
    @constraint(opt_model, eq, E * x .== d)
    @constraint(opt_model, ineq, M * x .<= h)
    @objective(opt_model, Min, 0.5 * x' * Q * x + c' * x + 0.5 * n)

    # apply the modifications, if any (from COBREXA, only optimizer based mods though)
    for mod in modifications
        mod(nothing, opt_model)
    end

    optimize!(opt_model)

    @assert(termination_status(opt_model) == JuMP.OPTIMAL)

    #: differentiate
    z = value.(opt_model[:x])
    n_vars = length(z)
    ν = dual.(opt_model[:eq])
    n_ν = length(ν)
    λ = dual.(opt_model[:ineq])
    n_λ = length(λ)

    vars = [z; ν; λ]

    if use_analytic
        n_mets = size(Ef(θ), 1)
        n_cons = size(M, 1)

        A = [
            Q -Ef(θ)' -M'
            Ef(θ) zeros(n_mets, n_ν) zeros(n_mets, n_λ)
            diagm(λ)*M zeros(n_cons, n_ν) diagm(M * z - hf(θ))
        ]
        B = dFdθ(z, ν, λ, θ)
    else
        #: KKT function
        F(x, θ) = Array(
            [
                Q * x[1:n_vars] + c - Ef(θ)' * x[n_vars.+(1:n_ν)] -
                M' * x[(n_vars+n_ν).+(1:n_λ)]
                Ef(θ) * x[1:n_vars] - d
                diagm(x[(n_vars+n_ν).+(1:n_λ)]) * (M * x[1:n_vars] - hf(θ))
            ],
        )

        A = ForwardDiff.jacobian(x -> F(x, θ), vars)
        B = ForwardDiff.jacobian(θ -> F(vars, θ), θ)
    end

    dx = -sparse(A) \ B #! will fail if det(A) = 0
    dx = dx[1:n_vars, :] # only return derivatives of variables

    # normalize flux direction
    if scale_output
        ndx = similar(dx)
        for i = 1:size(dx, 1)
            for j = 1:size(dx, 2)
                ndx[i, j] = round(θ[j] / vars[i] * dx[i, j]; digits = 8)
            end
        end
        return value.(x), ndx, objective_value(opt_model)
    else
        return value.(x), dx, objective_value(opt_model)
    end
end

"""
"""
function qp_objective_measured(
    rids,
    gids,
    obs_v_dict,
    obs_e_dict;
    vtol = 1e-3,
    etol = 1e-3,
    reg = 1e-1,
)
    n_vars = length(gids) + length(rids)
    c = zeros(n_vars)
    q = zeros(n_vars)
    n = 0

    for (i, rid) in enumerate(rids)
        if !haskey(obs_v_dict, rid) || abs(obs_v_dict[rid]) < vtol
            q[i] = reg
        else
            c[i] = -1.0 / abs(obs_v_dict[rid]) #! fluxes are positive in model
            q[i] = 1.0 / obs_v_dict[rid]^2
            n += 1
        end
    end
    k = length(rids)
    for (i, gid) in enumerate(gids)
        if !haskey(obs_e_dict, gid) || abs(obs_e_dict[gid]) < etol
            q[k+i] = reg
        else
            c[k+i] = -1.0 / obs_e_dict[gid]
            q[k+i] = 1.0 / obs_e_dict[gid]^2
            n += 1
        end
    end

    return spdiagm(q), sparse(c), n
end


