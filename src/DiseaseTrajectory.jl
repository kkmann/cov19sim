struct DiseaseTrajectory{T}
    vl::Vector{T}
    symptomatic::Vector{Bool}
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

duration(dt::DiseaseTrajectory) = length(dt.vl)
