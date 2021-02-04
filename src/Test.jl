abstract type Test end

length(x::Test) = 1
iterate(x::Test) = (x, nothing)
iterate(x::Test, state) = nothing

type(test::Test) = test.type

# limit of detection on VL scale
get_lod(test::Test) = throw(MethodError("not implemented"))

function get_probability_positive(test::Test, indv::Individual)
    if get_viral_load(indv) <= get_lod(test)
        1 - specificity(test, indv)
    else
        sensitivity(test, indv)
    end
end

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

get_lod(test::LogPropTest{T}) where {T<:Real} = 10 .^ test.log10_lod

function sensitivity(test::LogPropTest{T1}, indv::T2) where {T1<:Real,T2<:Individual}
    viral_load = get_viral_load(indv)
    elog10vl = max(0.0, log(10, viral_load) - test.log10_lod)
    min(1, max(0, test.slope * elog10vl))
end

specificity(test::LogPropTest{T1}, indv::T2) where {T1<:Real,T2<:Individual} = test.specificity



struct LogGLMTest{T} <: Test
    type::String
    lod::T
    beta_vl::T
    beta_ranef::T
    intercept::T
    specificity::T
end

get_lod(test::LogGLMTest{T}) where {T<:Real} = test.lod
specificity(test::LogGLMTest{T1}, indv::T2) where {T1<:Real,T2<:Individual} = test.specificity

function sensitivity(test::LogGLMTest{T1}, indv::T2) where {T1<:Real,T2<:Individual}
    elog10vl = log(10, get_viral_load(indv) - test.lod)
    if elog10vl <= 0
        return 1 - test.specificity
    else
        return max(1 - test.specificity, 1 / ( 1 + exp(-(
            test.beta_vl*elog10vl + test.beta_ranef*indv.dt.u + test.intercept
        ))))
    end
end
