import Pkg; Pkg.activate("..")

using Documenter, cov19sim

push!(LOAD_PATH,"../src/")
makedocs(sitename = "cov19sim")
