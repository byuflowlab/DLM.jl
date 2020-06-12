using Documenter, DLM

makedocs(;
    modules = [DLM],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Example" => "example.md",
        "Library" => "library.md"
    ],
    sitename = "DLM.jl",
    authors = "Taylor McDonnell <taylormcd@byu.edu>",
)

deploydocs(
    repo="github.com/byuflowlab/DLM.jl.git",
)
