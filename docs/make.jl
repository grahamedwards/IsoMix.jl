using IsoMix
using Documenter

DocMeta.setdocmeta!(IsoMix, :DocTestSetup, :(using IsoMix); recursive=true)

makedocs(;
    modules=[IsoMix],
    authors="Graham Harper Edwards",
    sitename="IsoMix.jl",
    format=Documenter.HTML(;
        canonical="https://grahamedwards.github.io/IsoMix.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/grahamedwards/IsoMix.jl",
    devbranch="main",
)
