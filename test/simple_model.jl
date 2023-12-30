
#=
Copyright (c) 2023, Heinrich-Heine University Duesseldorf
Copyright (c) 2023, University of Luxembourg

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Changes from copied code are indicated.
=#

import AbstractFBCModels as AM

# build simple model 
mets = Dict(
    "m1" => AM.CanonicalModel.Metabolite(),
    "m2" => AM.CanonicalModel.Metabolite(),
    "m3" => AM.CanonicalModel.Metabolite(),
    "m4" => AM.CanonicalModel.Metabolite(),
)

rxns = Dict(
    "r1" => AM.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m1" => 1.0),
    ),
    "r2" => AM.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 1.0,
        stoichiometry = Dict("m2" => 1.0),
    ),
    "r3" => AM.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m1" => -1.0, "m2" => -1.0, "m3" => 1.0),
    ),
    "r4" => AM.CanonicalModel.Reaction(
        lower_bound = -100.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m3" => -1.0, "m4" => 1.0),
    ),
    "r5" => AM.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m2" => -1.0, "m4" => 1.0),
    ),
    "r6" => AM.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m4" => -1.0),
        objective_coefficient = 1.0,
    ),
)

gs = Dict("g$i" => AM.CanonicalModel.Gene() for i = 1:5)

model = AM.CanonicalModel.Model(rxns, mets, gs)
