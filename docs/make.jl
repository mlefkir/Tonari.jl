push!(LOAD_PATH, "../src/")

using Documenter, Tonari
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"), style = :authoryear)

makedocs(sitename = "Tonari.jl",
	pages = ["Home" => "index.md",
		"Simulations" => ["simulations.md"],
        "Spectral analysis" => ["periodogram.md"],
		"References" => "references.md", 
        "API Reference" => "api.md",
	], format = Documenter.HTML(description = "Tonari.jl:", prettyurls = true), plugins = [bib])

deploydocs(
	repo = "github.com/mlefkir/Tonari.jl.git",
)
