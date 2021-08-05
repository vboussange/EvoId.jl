using Documenter, EvoId
# push!(LOAD_PATH,"/Users/victorboussange/ETHZ/projects/EvoId/") # not sure this is necessary
pathsrc = joinpath(@__DIR__,"src")
makedocs(sitename="EvoId.jl",
        format = Documenter.HTML(prettyurls = false),
        authors = "Victor Boussange",
        pages = [
            "Home" => "index.md",
            "Manual" => ["manual/agent.md", 
                        "manual/space.md", 
                        "manual/world.md", 
                        "manual/run_world.md"],
            "Examples" => [ joinpath(s[end-1:end]...) for s in splitpath.(readdir(joinpath(pathsrc,"examples"),join=true))],
            "Developping" => [ joinpath(s[end-1:end]...) for s in splitpath.(readdir(joinpath(pathsrc,"dev"),join=true))],
        # "contributing.md",
        ],)

deploydocs(repo = "github.com/vboussange/EvoId.jl")
