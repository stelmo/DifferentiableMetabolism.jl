# DifferentiableMetabolism.jl

[docs-img-stable]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url-stable]: https://stelmo.github.io/DifferentiableMetabolism.jl

[docs-img-dev]: https://img.shields.io/badge/docs-latest-0af.svg
[docs-url-dev]: https://stelmo.github.io/DifferentiableMetabolism.jl/dev

[ci-img]: https://github.com/stelmo/DifferentiableMetabolism.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/stelmo/DifferentiableMetabolism.jl/actions/workflows/ci.yml

[cov-img]: https://codecov.io/gh/stelmo/DifferentiableMetabolism.jl/branch/master/graph/badge.svg?token=A2ui7exGIH
[cov-url]: https://codecov.io/gh/stelmo/DifferentiableMetabolism.jl

[repostatus-url]: https://www.repostatus.org/#active
[repostatus-img]: https://www.repostatus.org/badges/latest/active.svg

| **Documentation** | **Tests** | **Coverage** | **Project status** |
|:---:|:---:|:---:|:---:|
| [![docs-img-stable]][docs-url-stable] [![docs-img-dev]][docs-url-dev] | [![CI status][ci-img]][ci-url] | [![codecov][cov-img]][cov-url] | [![repostatus-img]][repostatus-url] |


This package extends [COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl)
with the ability to differentiate an optimal solution of a constraint-based
metabolic model with respect to parameters.

To use this package, [download and install Julia](https://julialang.org/downloads/), and add
the following packages using the built in package manager:
```julia
] add COBREXA, DifferentiableMetabolism
```
Any optimization solver that is compatible with [JuMP](https://jump.dev/)
can be used, provided it can solve the optimization problems you are interested
in (typically, LPs and QPs), and it returns the dual variables through the JuMP interface.
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
