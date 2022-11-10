using RelativisticDynamics
using Documenter

DocMeta.setdocmeta!(RelativisticDynamics, :DocTestSetup, :(using RelativisticDynamics); recursive=true)

makedocs(;
    modules=[RelativisticDynamics],
    authors="Tom Kimpson",
    repo="https://github.com/tomkimpson/RelativisticDynamics.jl/blob/{commit}{path}#{line}",
    sitename="RelativisticDynamics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tomkimpson.github.io/RelativisticDynamics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "How to run RelativisticDynamics.jl"=>"how_to_run.md",
        "Initial Conditions"=>"IC.md",
        "MPD Equations"=>"MPD.md",
        "Visualisation"=>"visualisation.md",
        "Functions Index"=>"functions_index.md",
    ],
)

deploydocs(;
    repo="github.com/tomkimpson/RelativisticDynamics.jl",
    devbranch="main",
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"]
)
