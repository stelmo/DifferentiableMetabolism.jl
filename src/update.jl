"""
    update_Q!(diffmodel::DifferentiableModel, Q)

Update Q.
"""
update_Q!(diffmodel::DifferentiableModel, Q) = diffmodel.Q = Q

"""
    update_c!(diffmodel::DifferentiableModel, c)

Update c.
"""
update_c!(diffmodel::DifferentiableModel, c) = diffmodel.c = c

"""
    update_E!(diffmodel::DifferentiableModel, E)

Update E.
"""
update_E!(diffmodel::DifferentiableModel, E) = diffmodel.E = E

"""
    update_d!(diffmodel::DifferentiableModel, d)

Update d.
"""
update_d!(diffmodel::DifferentiableModel, d) = diffmodel.d = d

"""
    update_M!(diffmodel::DifferentiableModel, M)

Update M.
"""
update_M!(diffmodel::DifferentiableModel, M) = diffmodel.M = M

"""
    update_h!(diffmodel::DifferentiableModel, h)

Update h.
"""
update_h!(diffmodel::DifferentiableModel, h) = diffmodel.h = h
