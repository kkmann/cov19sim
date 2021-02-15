# taken from "Nosferican" https://discourse.julialang.org/t/function-like-expand-grid-in-r/4350/18
function gridexpand(arrays::AbstractVecOrMat...)
    l = size.(arrays, 1)
    nrows = prod(l)
    output = mapreduce(a_o -> repeat(a_o[1],
                                     inner = (a_o[2], 1),
                                     outer = (div(nrows, size(a_o[1], 1) * a_o[2]), 1)),
                       hcat,
                       zip(arrays, cumprod(prepend!(collect(l[1:end - 1]), 1))))
    return output
end

is_workday(day; workdays = collect(0:4)) = mod(day, 7) in workdays

inverse_logit(x::T) where {T<:Real} = 1 / (1 + exp(-x))
