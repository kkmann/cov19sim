abstract type Test end

length(x::Test) = 1
iterate(x::Test) = (x, nothing)
iterate(x::Test, state) = nothing

type(test::Test) = test.type

specificity(test::T) where {T<:Test} = test.specificity # default

sensitivity(test::T1, vl::T2) where {T1<:Test,T2<:Real} = throw(MethodError("not implemented"))
sensitivity(test::T, indv::I) where {T<:Test,I<:Individual} = sensitivity(test, get_viral_load(indv))

get_ar_window(test::T) where{T<:Test} = test.ar_window
get_ar_coefficient(test::T) where{T<:Test} = test.ar_coefficient

get_probability_positive(test::T, vl::R) where {T<:Test,R<:Real} = max( 1 - specificity(test), sensitivity(test, vl) )
function get_probability_positive(test::T, indv::I) where {T<:Test,I<:Individual}
    pr = get_probability_positive(test, get_viral_load(indv))
    n_log = length(indv.test_log_type)
    k = get_ar_window(test)
    if n_log > 0
        # look thorugh logs backwards till last test of same type
        for i in n_log:-1:1
            if now(indv) - indv.test_log_day[i] <= k
                if indv.test_log_type[i] == type(test)
                    a = get_ar_coefficient(test)
                    pr = (1 - a)*pr + a*indv.test_log_result[i]
                    break # found it, stop loop
                end
            else
                break # test log is sorted in time, stop searching backwards
            end
        end
    end
    pr
end

compliant(test::Test, indv::Individual) = type(test) == "pcr" ? true : compliant(indv)

function conduct_test!(test::Test, indv::Individual)
    pr = get_probability_positive(test, indv)
    if !compliant(test, indv)
        res = 2 # void
    else
        res = Int(rand(Distributions.Bernoulli(pr)))
    end
    log_test!(indv, type(test), res, pr)
    res
end
is_positive(result::Int)::Bool = result == 1
is_negative(result::Int)::Bool = result == 0
is_noncompliant(result::Int)::Bool = result == 3



struct FixedTest{T} <: Test
    type::String
    lod::T
    sensitivity::T
    specificity::T
    ar_window::Int
    ar_coefficient::T
end
FixedTest(type::String, sensitivity::T, specificity::T; lod::T = 0.0, ar_window::Int = 0, ar_coefficient::T = 0.0) where {T<:Real} = FixedTest{T}(type, lod, sensitivity, specificity, ar_window, ar_coefficient)
standard_pcr_test = FixedTest("pcr", .975, 1.0; lod = 150.0)

sensitivity(test::FixedTest{T1}, vl::T2) where {T1<:Real,T2<:Real} = vl >= test.lod ? test.sensitivity : 0.0



struct LogRegTest{T} <: Test
    type::String
    slope::T
    intercept::T
    specificity::T
    ar_window::Int
    ar_coefficient::T
end
LogRegTest(type::String, slope::T, intercept::T, specificity::T; ar_window::Int = 0, ar_coefficient::T = 0.0) where {T<:Real} =
    LogRegTest{T}(type, slope, intercept, specificity, ar_window, ar_coefficient)

function sensitivity(test::LogRegTest{T1}, vl::T2) where {T1<:Real,T2<:Real}
    inverse_logit( test.slope*log10(vl) + test.intercept )
end
