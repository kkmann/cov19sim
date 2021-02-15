abstract type InfectionModel end

function get_infection_probability(im::InfectionModel, dt::DiseaseTrajectory{T1}, day::T2)::T1 where
    {T1<:Real,T2<:Int}

    throw(MethodError("not implemented"))
end



struct LogPropInfectionModel{T} <: InfectionModel
    gamma::T
    lli::T
end

function get_infection_probability(im::LogPropInfectionModel{T1}, dt::DiseaseTrajectory{T1}, day::T2)::T1 where
    {T1<:Real,T2<:Int}

    vl = get_viral_load(dt, day)
    if vl <= lli
        return 0.0
    end
    min(1, max(0, im.gamma * log(10, vl) ) )
end



struct LogRegInfectionModel{T} <: InfectionModel
    gamma::T
    lli::T
    pr_lli::T
end
LogRegInfectionModel(gamma::T, lli::T; pr_lli::T = 0.001) where {T<:Real} = LogRegInfectionModel{T}(gamma, lli, pr_lli)

function get_infection_probability(im::LogRegInfectionModel{T1}, dt::DiseaseTrajectory{T1}, day::T2)::T1 where
    {T1<:Real,T2<:Int}

    vl = get_viral_load(dt, day)
    if vl <= 0
        return 0.0
    end
    intercept = log( im.pr_lli / (1 - im.pr_lli) )
    linpred = im.gamma * log(10, max(0, vl - im.lli)) # linear predictor of logreg
    inverse_logit(linpred + intercept)
end
