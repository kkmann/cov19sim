# abstract type declared in Group.jl

function test_and_isolate!(pol::T1, gr::Group) where {T1<:Policy}
    throw(MethodError("not implemented"))
end

length(x::Policy) = 1
iterate(x::Policy) = (x, nothing)
iterate(x::Policy, state) = nothing
