using Documenter, ABMEv
# push!(LOAD_PATH,"/Users/victorboussange/ETHZ/projects/ABMEv/") # not sure this is necessary
pathsrc = joinpath(@__DIR__,"src")
makedocs(sitename="ABMEv.jl",
        format = Documenter.HTML(prettyurls = false),
        pages = [
            "Home" => "index.md",
            "Manual" => [ joinpath(s[end-1:end]...) for s in splitpath.(readdir(joinpath(pathsrc,"manual"),join=true))],
            "Mathematics" => [ joinpath(s[end-1:end]...) for s in splitpath.(readdir(joinpath(pathsrc,"mathematics"),join=true))],
            "Examples" => [ joinpath(s[end-1:end]...) for s in splitpath.(readdir(joinpath(pathsrc,"examples"),join=true))],
            "Library" => Any[
                "Public" => "lib/public.md",
                # "Internals" => map(
                #     s -> "lib/internals/$(s)",
                #     sort(readdir(joinpath(@__DIR__, "src/lib/internals")))
                # ),
                ],
            "Developping" => [ joinpath(s[end-1:end]...) for s in splitpath.(readdir(joinpath(pathsrc,"dev"),join=true))],
        # "contributing.md",
        ],)

deploydocs(repo = "github.com/vboussange/ABMEv.jl",
            deps = Deps.pip("mkdocs", "python-markdown-math"),)
