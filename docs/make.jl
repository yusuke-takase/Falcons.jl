using Falcons
using Documenter

DocMeta.setdocmeta!(Falcons, :DocTestSetup, :(using Falcons); recursive = true)

makedocs(;
	modules = [Falcons],
	authors = "Yusuke Takase",
	repo = Remotes.GitHub("yusuke-takase", "Falcons.jl"),
	sitename = "Falcons.jl",
	format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true", assets = String[]),
	pages = [
		"Home" => "index.md",
		"Conventions" => "conventions.md",
		"Time-dependent rotations" => "attitude.md",
		"Instrument database" => "imo.md",
		"GPU backend" => "gpu.md",
		"Migrating from 1.0" => "migration.md",
		"API reference" => "api.md",
	],
	checkdocs = :exports,
)
