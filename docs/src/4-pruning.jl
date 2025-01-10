
import DifferentiableMetabolism as D
import FastDifferentiation as F
const Ex = F.Node
import ConstraintTrees as C
import AbstractFBCModels as A
import JSONFBCModels as JFBC
import COBREXA as X
import Tulip as T

cd("docs")
include("../../test/data_static.jl")

# Load model, and convert to CanonicalModel for ease of use
model = convert(A.CanonicalModel.Model, X.load_model("e_coli_core.json"))

for rid in ["ACt2r", "ETOHt2r", "PYRt2", "SUCCt3"]
    model.reactions[rid].gene_association_dnf = [["s0001"]]
end

# Modify the model a little bit
model.reactions["EX_glc__D_e"].lower_bound = -1000.0 # capacity bound suffices
model.reactions["PFL"].upper_bound = 0.0 # aerobic simulation

reaction_isozymes = Dict{String,Dict{String,X.Isozyme}}()

# Populate reaction_isozymes with parameters
for rid in A.reactions(model)
    grrs = A.reaction_gene_association_dnf(model, rid)
    isnothing(grrs) && continue # skip if no grr available
    haskey(ecoli_core_reaction_kcats, rid) || continue # skip if no kcat data available

    for (i, grr) in enumerate(grrs)

        kcat = ecoli_core_reaction_kcats[rid] * 3.6 # change unit to k/h

        d = get!(reaction_isozymes, rid, Dict{String,X.Isozyme}()) 
        d["isozyme_$i"] = X.Isozyme( 
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), 
            kcat_forward = kcat, 
            kcat_reverse = kcat, 
        ) 
    end
end

gene_product_molar_masses = Dict(k => v for (k, v) in ecoli_core_gene_product_masses)


m = X.enzyme_constrained_flux_balance_constraints(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = 50.0,
)

sol = X.optimized_values(
    m,
    objective = m.objective.value,
    optimizer = T.Optimizer,
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

vs = D.variable_order(m)
zero_tol = 1e-6

zero_idxs = Set{Int}() # all variables that are close to zero

C.itraverse(m) do p, x
    if length(x.value.idxs) == 1 && abs(foldl(getproperty, p; init=sol)) <= zero_tol # TODO should have gene and flux specific zeros?
        push!(zero_idxs, first(x.value.idxs))
    end
end

vars = [idx in zero_idxs ? zero(C.LinearValue) : C.variable(; idx).value for idx = 1:C.variable_count(m)]

pm = C.substitute(m, vars)

pm = C.filter_leaves(pm) do leaf
    !isempty(leaf.value.idxs)
end

C.variable_count(pm) # TODO delete

# now get rid of forward/reverse fluxes
C.itraverse(m) do p, x
    if first(p) == :fluxes
        if foldl(getproperty, p; init=sol) > 0
            x.bound = C.Between(0,1000)
        else
            x.bound = C.Between(-1000,0)
        end
    end
end





pm = C.drop_zeros(pm)
pm = C.prune_variables() # CTs is wonderful


C.variable_count(pm) # TODO delete

parameters = Set{Symbol}()
C.traverse(pm) do leaf
    ps = []
    leaf.value isa LinearValueP && append!(ps, leaf.value.weights)
    leaf.value isa QuadraticValueP && append!(ps, leaf.value.weights)
    leaf.bound isa EqualToP && push!(ps, leaf.bound.equal_to)
    leaf.bound isa BetweenP && append!(ps, [leaf.bound.lower, leaf.bound.upper])

    for s in filter(!isreal, F.value.(ps))
        push!(parameters, s)
    end
end
