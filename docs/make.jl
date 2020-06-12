using Documenter
using DLM

makedocs(
    sitename = "DLM",
    pages = [
        "Home" => "index.md",
        "Example" => "example.md",
        "Library" => "library.md"
    ],
#    format = Documenter.HTML(),
    modules = [DLM]
)

# deploydocs(
#     repo="github.com/byuflowlab/DLM.jl.git",
# )
