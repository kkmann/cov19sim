abstract type Test end

length(x::Test) = 1
iterate(x::Test) = (x, nothing)
iterate(x::Test, state) = nothing

type(test::Test) = test.type

specificity(test::T) where {T<:Test} = test.specificity # default

sensitivity(test::T1, vl::T2, u::T3) where {T1<:Test,T2<:Real,T3<:Real} = throw(MethodError("not implemented"))
sensitivity(test::T, indv::I) where {T<:Test,I<:Individual} = sensitivity(test, get_viral_load(indv), indv.dt.u)

function get_probability_positive(test::T, vl::F1, u::F2) where {T<:Test,F1<:Real,F2<:Real}
   max(
       1 - specificity(test),
       sensitivity(test, vl, u)
    )
end
get_probability_positive(test::T, indv::I) where {T<:Test,I<:Individual} = get_probability_positive(test, get_viral_load(indv), indv.dt.u)

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

sensitivity(test::FixedTest{T1}, vl::T2, u::T3) where {T1<:Real,T2<:Real,T3<:Real} = test.sensitivity



struct LogRegTest{T} <: Test
    type::String
    beta_vl::T
    beta_u::T
    intercept::T
    specificity::T
end

function sensitivity(test::LogRegTest{T1}, vl::T2, u::T3) where {T1<:Real,T2<:Real,T3<:Real}
    inverse_logit( test.beta_vl*log(10, vl) + test.beta_u*u + test.intercept )
end
