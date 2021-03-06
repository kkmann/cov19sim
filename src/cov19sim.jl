module cov19sim

import Distributions
import Interpolations
import DataFrames, DataFramesMeta, Statistics

import LinearAlgebra.Symmetric

import Base.length, Base.iterate, Base.UUID
import UUIDs.uuid4
length(x::UUID) = 1
iterate(x::UUID) = (x, nothing)
iterate(x::UUID, state) = nothing

include("util.jl")

include("DiseaseTrajectory.jl")
export DiseaseTrajectory
export get_viral_load, is_symptomatic, has_recovered, duration

include("InfectionModel.jl")
export InfectionModel, ProportionalInfectionModel, LogRegInfectionModel, LogPropInfectionModel
export get_infection_probability

include("DiseaseModel.jl")
export DiseaseModel
export sample
include("disease-models/LarremoreModel.jl")
export LarremoreModel
include("disease-models/HeavyTailsModel.jl")
export HeavyTailsModel

include("Individual.jl")
export Individual
export now, step!, steps!, infect!, meet!, get_status_log, get_contact_log, log_test!, get_test_log,
    get_status_logs, get_contact_logs, get_test_logs, is_infected,
    is_pcr_positive, is_newly_symptomatic

include("Test.jl")
export Test
export conduct_test!, type, sensitivity, specificity, get_probability_positive
export FixedTest, LogRegTest

include("Group.jl")
export Group, Policy
export test_and_isolate!

include("Policy.jl")
include("policies/DoNothing.jl")
export DoNothing
include("policies/SymptomaticIsolation.jl")
export  SymptomaticIsolation
include("policies/DynamicScreening.jl")
export DynamicScreening

include("Population.jl")
export Population
export get_status_over_time, get_test_logs_over_time, get_mean_contacts_over_time, evaluate,
n_individuals, n_infected, n_workdays_missed, n_tests, n_infectious_per_day,
mean_contacts_per_day, get_adjacency_matrix, resample

include("ThreeLevelPopulation.jl")
export ThreeLevelPopulation

end # module
