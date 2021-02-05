struct DiseaseTrajectory{T}
    vl::Vector{T}
    symptomatic::Vector{Bool}
    u::T # random effect for test sensitivity
end

length(dt::DiseaseTrajectory) = 1
iterate(dt::DiseaseTrajectory) = (dt, nothing)
iterate(dt::DiseaseTrajectory, state) = nothing

function get_viral_load(DT::DiseaseTrajectory{T1}, day::T2) where {T1<:Real,T2<:Int}
    (day < 0) | (day >= length(DT.vl)) ? 0.0 : DT.vl[day + 1]
end

function is_symptomatic(DT::DiseaseTrajectory{T1}, day::T2) where {T1<:Real,T2<:Int}
    (day < 0) | (day >= length(DT.vl)) ? false : DT.symptomatic[day + 1]
end

function has_recovered(DT::DiseaseTrajectory{T1}, day::T2) where {T1<:Real,T2<:Int}
    day >= length(DT.vl)
end

"""
    duration(dt::DiseaseTrajectory)

Returns the duration (in days) of the disease trajecotry `dt`, the duration is the number of days
from the infection until the final clearance (viral load zero).
"""
duration(dt::DiseaseTrajectory) = length(dt.vl)
