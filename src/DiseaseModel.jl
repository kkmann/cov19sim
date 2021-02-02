abstract type DiseaseModel end

length(dm::DiseaseModel) = 1
iterate(dm::DiseaseModel) = (dm, nothing)
iterate(dm::DiseaseModel, state) = nothing

function (get_infection_probability(DM::T1, DT::DiseaseTrajectory{T2}, day::T3)::T2) where
    {T1<:DiseaseTrajectory,T2<:Real,T3<:Int}

    throw(MethodError("not implemented"))
end

function (sample(DM::T1)::DiseaseTrajectory{T2}) where {T1<:DiseaseTrajectory,T2<:Real}

    throw(MethodError("not implemented"))
end
