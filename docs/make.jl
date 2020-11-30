using Documenter, BDisposal

push!(LOAD_PATH,"../src/")
makedocs(sitename="BDisposal.jl Documentation",
         #root = "../",
         #source = "src",
         #build = "build",
         pages = [
            "Index" => "index.md",
            "Efficiency Scores" => "efficiencyScores.md",
            "Productivity Index" => "prodIndex.md",
         ],
         format = Documenter.HTML(prettyurls = false)
)
deploydocs(
    repo = "github.com/sylvaticus/BDisposal.jl.git",
    devbranch = `main`
)
