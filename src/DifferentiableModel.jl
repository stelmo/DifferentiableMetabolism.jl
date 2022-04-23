"""
    mutable struct DifferentiableModel

A struct representing a generalized differentiable constraint based model.
Format of the model is:
```
minimize    ½ * x' * Q(θ) * x + c(θ)' * x 
s.t.        E(θ) * x = d(θ) 
            M(θ) * x ≤ h(θ)
```
Nominally, every component can be differentiated. Analytic derivatives of the
KKT function with respect to the variables and parameters can be used to speed
up the derivative calculation. For this, the derivative of `θ` with respect to
the KKT conditions should be supplied through `param_derivs`, which should be 
a function taking 4 arguments, `(x, ν, λ, θ)`.

# Fields
```
Q::Function
c::Function
E::Function
d::Function
M::Function
h::Function
θ::Vector{Float64}
auto_derivs::Function
param_derivs::Function
param_ids::Vector{String}
var_ids::Vector{String}
```
"""
mutable struct DifferentiableModel
    Q::Function
    c::Function
    E::Function
    d::Function
    M::Function
    h::Function
    θ::Vector{Float64}
    auto_derivs::Function
    param_derivs::Function
    param_ids::Vector{String}
    var_ids::Vector{String}
end

"""
Pretty printing of a differentiable model.
"""
function Base.show(io::IO, ::MIME"text/plain", m::DifferentiableModel)
    println(io, "", "Differentiable metabolic with: ")
    println(io, "    ", length(m.var_ids), " variables")
    println(io, "    ", length(m.param_ids), " parameters")
end

