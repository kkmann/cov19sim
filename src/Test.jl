abstract type Test end

length(x::Test) = 1
iterate(x::Test) = (x, nothing)
iterate(x::Test, state) = nothing

type(test::Test) = test.type

get_probability_positive(test::Test, indv::Individual) = max(1 - specificity(test, indv), sensitivity(test, indv))

function conduct_test!(test::Test, indv::Individual)
    pr = get_probability_positive(test, indv)
    res = rand(Distributions.Bernoulli(pr))
    log_test!(indv, type(test), res, pr)
    res
end




struct LogPropTest{T} <: Test
    type::String
    log10_lod::T
    slope::T
    specificity::T
end

function sensitivity(test::LogPropTest{T1}, indv::T2) where {T1<:Real,T2<:Individual}
    viral_load = get_viral_load(indv)
    elog10vl = max(0.0, log(10, viral_load) - test.log10_lod)
    min(1, max(0, test.slope * elog10vl))
end

specificity(test::LogPropTest{T1}, indv::T2) where {T1<:Real,T2<:Individual} = test.specificity



struct LogRegTest{T} <: Test
    type::String
    beta_vl::T
    beta_u::T
    intercept::T
    specificity::T
end

specificity(test::LogRegTest{T1}, indv::T2) where {T1<:Real,T2<:Individual} = test.specificity

function sensitivity(test::LogRegTest{T1}, indv::T2) where {T1<:Real,T2<:Individual}
    vl = get_viral_load(indv)
    if vl <= 0
        return 1 - test.specificity
    end
    elog10vl = log(10, get_viral_load(indv))
    inverse_logit( test.beta_vl*elog10vl + test.beta_u*indv.dt.u + test.intercept )
end
