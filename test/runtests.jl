using cov19sim, Test, DataFrames, DataFramesMeta

import Random.seed!, Random.randperm
seed!(42)

#include("tests/test_individual.jl")
include("tests/test_tests.jl")
include("tests/test_Group.jl")
include("tests/test_population.jl")
include("tests/test_SymptomaticIsolation.jl")
include("tests/test_RegularScreening.jl")
include("tests/test_DynamicScreening.jl")
