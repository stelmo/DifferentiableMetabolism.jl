using Documenter, Literate, DifferentiableMetabolism

examples = sort(filter(x -> endswith(x, ".jl"), readdir(joinpath(@__DIR__, "src"), join = true)))

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
            assets=String[],
        ),
        sitename = "DifferentiableMetabolism.jl",
        linkcheck = false,
        pages = ["README" => "indeCOBREXA.md"; example_mds; "Reference" => "reference.md"],
    )
end

deploydocs(
    repo = "github.com/stelmo/DifferentiableMetabolism.jl",
    target = "build",
    branch = "gh-pages",
    push_preview = false,
)

