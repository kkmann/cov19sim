module cov19sim

import Distributions
import Interpolations
import DataFrames, DataFramesMeta, Statistics

import Base.length, Base.iterate, Base.UUID
import UUIDs.uuid4
length(x::UUID) = 1
iterate(x::UUID) = (x, nothing)
iterate(x::UUID, state) = nothing

include("DiseaseTrajectory.jl")
export DiseaseTrajectory
export get_viral_load, is_symptomatic, has_recovered

include("DiseaseModel.jl")
export DiseaseModel
export sample, get_infection_probability
include("LarremoreModel.jl")
export LarremoreModel

include("Individual.jl")
export Individual
export now, pcr_test_and_isolate!, step!, steps!, infect!, meet!, get_status_log, get_contact_log, log_test!, get_test_log,
    get_status_logs, get_contact_logs, get_test_logs, is_infected

include("Test.jl")
export Test
export conduct_test!, type, sensitivity, specificity
export LogPropTest

include("Group.jl")
export Group
export test_and_isolate!

include("policies.jl")
export DoNothing, SymptomaticIsolation, RegularScreening, DynamicScreening

include("Population.jl")
export Population
export get_status_over_time, evaluate

include("ThreeLevelPopulation.jl")
export ThreeLevelPopulation

end # module
