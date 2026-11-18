using Documenter
using BOLDTrace

DocMeta.setdocmeta!(BOLDTrace, :DocTestSetup, :(using BOLDTrace); recursive=true)

makedocs(;
    modules=[BOLDTrace],
    sitename="BOLDTrace.jl",
    authors="Alex Laini",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alexology.github.io/BOLDTrace.jl",
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "BOLD data and the LCA" => "lca.md",
        "Attribution, counters and citations" => "attribution.md",
        "Export files" => "exports.md",
        "API" => "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(; repo="github.com/alexology/BOLDTrace.jl", devbranch="main")
