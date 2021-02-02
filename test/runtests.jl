using cov19sim, Test

import Random.seed!

seed!(42)

#include("tests/test_individual.jl")
#include("tests/test_LogPropTest.jl")
# include("tests/test_Group.jl")
include("tests/test_ThreeLevelPopulation.jl")
