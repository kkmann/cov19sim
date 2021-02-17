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



struct FixedTest{T} <: Test
    type::String
    lod::T
    sensitivity::T
    specificity::T
end
FixedTest(type::String, sensitivity::T, specificity::T; lod::T = 0.0) where {T<:Real} = FixedTest{T}(type, lod, sensitivity, specificity)

standard_pcr_test = FixedTest("pcr", .975, 1.0; lod = 150.0)

function sensitivity(test::FixedTest{T1}, indv::T2) where {T1<:Real,T2<:Individual}
    if get_viral_load(indv) > test.lod
        return test.sensitivity
    else
        return 1 - test.specificity
    end
end
specificity(test::FixedTest{T1}, indv::T2) where {T1<:Real,T2<:Individual} = test.specificity



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
