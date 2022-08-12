using Falcons
using Documenter

DocMeta.setdocmeta!(Falcons, :DocTestSetup, :(using Falcons); recursive=true)

makedocs(;
    modules=[Falcons],
    authors="Yusuke Takase",
    repo="https://github.com/yusuke-takase/Falcons.jl/blob/{commit}{path}#{line}",
    sitename="Falcons.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yusuke-takase.github.io/Falcons.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Scanning" => "Scanning.md",
        #"Mapmake" => "Mapmake.md",
    ],
)

deploydocs(;
    repo="github.com/yusuke-takase/Falcons.jl",
)
