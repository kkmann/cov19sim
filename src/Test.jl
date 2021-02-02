abstract type Test end

length(x::Test) = 1
iterate(x::Test) = (x, nothing)
iterate(x::Test, state) = nothing

function conduct_test!(test::Test, indv::Individual)
    vl = get_viral_load(indv)
    if vl == 0.0
        pr = 1 - specificity(test)
    else
        pr = sensitivity(test, vl)
    end
    res = rand(Distributions.Bernoulli(pr))
    log_test!(indv, type(test), res)
    res
end

type(test::Test) = test.type




struct LogPropTest{T} <: Test
    type::String
    log10_lod::T
    slope::T
    specificity::T
end

function sensitivity(test::LogPropTest{T1}, viral_load::T2) where {T1<:Real,T2<:Real}
    elog10vl = max(0.0, log(10, viral_load) - test.log10_lod)
    min(1, max(0, test.slope * elog10vl))
end

specificity(test::LogPropTest{T1}) where {T1<:Real,T2<:Real} = test.specificity
