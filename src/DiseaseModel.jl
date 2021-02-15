abstract type DiseaseModel end

length(dm::DiseaseModel) = 1
iterate(dm::DiseaseModel) = (dm, nothing)
iterate(dm::DiseaseModel, state) = nothing

function (get_infection_probability(dm::T1, dt::DiseaseTrajectory{T2}, day::T3)::T2) where
    {T1<:DiseaseModel,T2<:Real,T3<:Int}

    get_infection_probability(dm.infection_model, dt, day)
end

function (sample(dm::T1)::DiseaseTrajectory{T2}) where {T1<:DiseaseTrajectory,T2<:Real}
    throw(MethodError("not implemented"))
end
