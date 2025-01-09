
#=
ExCopyright (c) 2025, Heinrich-Heine University Duesseldorf
ExCopyright (c) 2025, University of Luxembourg

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

makedocs(
    modules = [DifferentiableMetabolism],
    sitename = "DifferentiableMetabolism.jl",
    linkcheck = false,
    format = Documenter.HTML(
        ansicolor = true,
        canonical = "https://stelmo.github.io/DifferentiableMetabolism.jl/stable/",
    ),
    pages = [
        "README" => "index.md",
        "Examples" => [
            "Contents" => "examples.md"
            "examples_mds"
        ],
        "Reference" => "reference.md",
    ],
)

# clean up examples -- we do not need to deploy all the stuff that was
# generated in the process
#
# extra fun: failing programs (such as plotting libraries) may generate core
# dumps that contain the dumped environment strings, which in turn contain
# github auth tokens. These certainly need to be avoided.
examples_names = [n[begin:end-3] for n in examples]
ipynb_names = examples_names .* ".ipynb"
examples_allowed_files = vcat("index.html", ipynb_names)
@info "allowed files:" examples_allowed_files
for (root, dirs, files) in walkdir(joinpath(@__DIR__, "build", "examples"))
    for f in files
        if !(f in examples_allowed_files)
            @info "removing notebook build artifact `$(joinpath(root, f))'"
            rm(joinpath(root, f))
        end
    end
end

deploydocs(
    repo = "github.com/stelmo/DifferentiableMetabolism.jl",
    target = "build",
    branch = "gh-pages",
    devbranch = "master",
)
