"""
$(TYPEDSIGNATURES)

Update Q.
"""
update_Q!(diffmodel::DifferentiableModel, Q) = diffmodel.Q = Q

"""
$(TYPEDSIGNATURES)

Update c.
"""
update_c!(diffmodel::DifferentiableModel, c) = diffmodel.c = c

"""
$(TYPEDSIGNATURES)

Update E.
"""
update_E!(diffmodel::DifferentiableModel, E) = diffmodel.E = E

"""
$(TYPEDSIGNATURES)

Update d.
"""
update_d!(diffmodel::DifferentiableModel, d) = diffmodel.d = d

"""
$(TYPEDSIGNATURES)

Update M.
"""
update_M!(diffmodel::DifferentiableModel, M) = diffmodel.M = M

"""
$(TYPEDSIGNATURES)

Update h.
"""
update_h!(diffmodel::DifferentiableModel, h) = diffmodel.h = h
