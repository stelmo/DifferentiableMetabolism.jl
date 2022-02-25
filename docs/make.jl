using DifferentiableMetabolism
using Documenter

DocMeta.setdocmeta!(
    DifferentiableMetabolism,
    :DocTestSetup,
    :(using DifferentiableMetabolism);
    recursive = true,
)

makedocs(;
    modules = [DifferentiableMetabolism],
    authors = "St. Elmo Wilken <stelmozors@gmail.com> and contributors",
    repo = "https://github.com/stelmo/DifferentiableMetabolism.jl/blob/{commit}{path}#{line}",
    sitename = "DifferentiableMetabolism.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://stelmo.github.io/DifferentiableMetabolism.jl",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/stelmo/DifferentiableMetabolism.jl", devbranch = "master")
