"""
    mutable struct DifferentiableModel

A struct representing a generalized differentiable constraint based model.
Format of the model is:
```
minimize    ½ * x' * Q(θ) * x + c(θ)' * x 
s.t.        E(θ) * x = d(θ) 
            M(θ) * x ≤ h(θ)
```
Nominally, every component can be differentiated.

# Fields
```
Q::Function
c::Function
E::Function
d::Function
M::Function
h::Function
θ::Vector{Float64}
analytic_var_derivs::Function
analytic_par_derivs::Function
var_ids::Vector{String}
param_ids::Vector{String}
```
Note, to ensure differentiability, preprocessing of the model used to derive the
[`DifferentiableModel`](@ref) is required. In short, only an active solution may
be differentiated, this requires that:
- the base model does not possess any isozymes, each reaction may be catalyzed
  by one enzyme (complex) only,
- all the reactions should be unidirectinal,
- the kcats in `rid_enzyme` are for the appropriate direction used in the model,
- all rids in `rid_enzyme` are used in the model.
"""
mutable struct DifferentiableModel
    Q::Function
    c::Function
    E::Function
    d::Function
    M::Function
    h::Function
    θ::Vector{Float64}
    analytic_var_derivs::Function
    analytic_par_derivs::Function
    var_ids::Vector{String}
    param_ids::Vector{String}
end

"""
Pretty printing of a differentiable model.
"""
function Base.show(io::IO, ::MIME"text/plain", m::DifferentiableModel)
    println(io, "", "Differentiable metabolic with: ")
    println(io, "    ", length(m.var_ids), " variables")
    println(io, "    ", length(m.param_ids), " parameters")
end

