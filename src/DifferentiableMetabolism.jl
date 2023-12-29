module DifferentiableMetabolism

using DocStringExtensions

import COBREXA as X
import JuMP as J
import ConstraintTrees as C
import Symbolics as S
import LinearAlgebra as L

# extend constraint trees to make things differentiable
include("parameter_bound.jl")
include("parameter_linearvalue.jl")
include("parameter_isozyme.jl")
include("parameter_misc.jl")

end
