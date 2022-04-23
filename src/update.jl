function update_auto_derivs(diffmodel::DifferentiableModel)
    θ = diffmodel.θ
    diffmodel.auto_derivs = get_auto_derivs(
        diffmodel.Q,
        diffmodel.c,
        diffmodel.E,
        diffmodel.d,
        diffmodel.M,
        diffmodel.h,
        length(diffmodel.c(θ)),
        size(diffmodel.E(θ), 1),
        size(diffmodel.M(θ), 1),
        length(θ),
    )
    nothing
end

"""
    update_Q!(diffmodel::DifferentiableModel, Q)

Update Q and the automatic derivatives.
"""
function update_Q!(diffmodel::DifferentiableModel, Q)
    # update Q
    diffmodel.Q = Q

    # update automatic derivative function   
    update_auto_derivs(diffmodel)

    nothing
end

"""
    update_c!(diffmodel::DifferentiableModel, c)

Update c and the automatic derivatives.
"""
function update_c!(diffmodel::DifferentiableModel, c)
    # update c
    diffmodel.c = c

    # update automatic derivative function
    update_auto_derivs(diffmodel)

    nothing
end

"""
    update_E!(diffmodel::DifferentiableModel, E)

Update E and the automatic derivatives.
"""
function update_E!(diffmodel::DifferentiableModel, E)
    # update E
    diffmodel.E = E

    # update automatic derivative function
    update_auto_derivs(diffmodel)

    nothing
end

"""
    update_d!(diffmodel::DifferentiableModel, d)

Update d and the automatic derivatives.
"""
function update_d!(diffmodel::DifferentiableModel, d)
    # update d
    diffmodel.d = d

    # update automatic derivative function
    update_auto_derivs(diffmodel)

    nothing
end

"""
    update_M!(diffmodel::DifferentiableModel, M)

Update M and the automatic derivatives.
"""
function update_M!(diffmodel::DifferentiableModel, M)
    # update M
    diffmodel.M = M

    # update automatic derivative function
    update_auto_derivs(diffmodel)

    nothing
end

"""
    update_h!(diffmodel::DifferentiableModel, h)

Update h and the automatic derivatives.
"""
function update_h!(diffmodel::DifferentiableModel, h)
    # update h
    diffmodel.h = h

    # update automatic derivative function
    update_auto_derivs(diffmodel)

    nothing
end
