# DifferentiableMetabolism.jl

| Build status | Documentation |
|:---:|:---:|
| ![CI status](https://github.com/stelmo/DifferentiableMetabolism.jl/workflows/CI/badge.svg) [![codecov](https://codecov.io/gh/stelmo/DifferentiableMetabolism.jl/branch/master/graph/badge.svg?token=A2ui7exGIH&style=for-the-badge)](https://codecov.io/gh/stelmo/DifferentiableMetabolism.jl) | [![Documentation](https://img.shields.io/badge/documentation-8e44ad?style=for-the-badge)](https://stelmo.github.io/DifferentiableMetabolism.jl/dev) |

This package extends [COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl)
with the ability to differentiate an optimal solution of a constraint-based
metabolic model with respect to parameters. 

Note, only non-degenerate (unique) solutions can be differentiated for the
derivatives to have a concrete interpretation. You need to ensure that your
solution is non-degenerate, otherwise you will only compute sub-gradients.
Builtin functionality (pruning) can help with this.

To use this package, [download and install Julia](https://julialang.org/downloads/), and add 
the following packages using the built in package manager:
```julia
] add COBREXA, DifferentiableMetabolism
```
Any optimization solver that is compatible with [JuMP](https://jump.dev/)
can be used, provided it can solve the optimization problems you are interested
in (LPs and QPs), and it returns the dual variables through the JuMP interface.
In the tests we use [Tulip.jl](https://github.com/ds4dm/Tulip.jl) for LPs, and
[Clarabel.jl](https://github.com/oxfordcontrol/Clarabel.jl) for QPs. Other
solvers, like Gurobi, work well, but they require a license (usually free for
academic use).

You can test the installation through:
```julia
] test DifferentiableMetabolism
```

For more information, please see the documentation.

This package is maintained and open for extensions. Feel free to discuss changes
and ideas via issues and pull requests.

#### Acknowledgements

`DifferentiableMetabolism.jl` was developed at Institute for Quantitative and
Theoretical Biology at Heinrich Heine University Düsseldorf
([qtb.hhu.de](https://www.qtb.hhu.de/en/)), and at the Luxembourg Centre for
Systems Biomedicine of the University of Luxembourg
([uni.lu/lcsb](https://www.uni.lu/lcsb)).

<img src="docs/src/assets/hhu.svg" alt="HHU logo" height="64px" style="height:64px; width:auto">   <img src="docs/src/assets/qtb.svg" alt="QTB logo" height="64px" style="height:64px; width:auto">   <img src="docs/src/assets/unilu.svg" alt="Uni.lu logo" height="64px">   <img src="docs/src/assets/lcsb.svg" alt="LCSB logo" height="64px">
