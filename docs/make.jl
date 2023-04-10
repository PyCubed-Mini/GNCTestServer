# Inside make.jl
push!(LOAD_PATH, "../src/")
using SatellitePlayground
using Documenter
makedocs(
    sitename="SatellitePlayground.jl",
    modules=[SatellitePlayground],
    pages=[
        "Home" => "index.md"
    ])
deploydocs(;
    repo="https://github.com/PyGNC/SatellitePlayground.jl"
)