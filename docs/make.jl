
#=
Copyright (c) 2024, Heinrich-Heine University Duesseldorf
Copyright (c) 2024, University of Luxembourg

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=#

using Documenter, Literate, DifferentiableMetabolism

examples =
    sort(filter(x -> endswith(x, ".jl"), readdir(joinpath(@__DIR__, "src"), join = true)))

for example in examples
    Literate.markdown(
        example,
        joinpath(@__DIR__, "src"),
        repo_root_url = "https://github.com/stelmo/DifferentiableMetabolism.jl/blob/master",
    )
end

example_mds = first.(splitext.(basename.(examples))) .* ".md"

withenv("COLUMNS" => 150) do
    makedocs(
        modules = [DifferentiableMetabolism],
        clean = false,
        format = Documenter.HTML(
            ansicolor = true,
            canonical = "https://stelmo.github.io/DifferentiableMetabolism.jl",
            assets = String[],
        ),
        sitename = "DifferentiableMetabolism.jl",
        linkcheck = false,
        pages = ["README" => "index.md"; example_mds; "Reference" => "reference.md"],
    )
end

deploydocs(
    repo = "github.com/stelmo/DifferentiableMetabolism.jl",
    target = "build",
    branch = "gh-pages",
    push_preview = false,
)
