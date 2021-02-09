
abstract type Policy end

mutable struct Group{T1,T2,T3}
    id::String
    individuals::Vector{T1}
    policy::T2
    meeting_probability::T3
    meeting_days::Vector{Int}
end

length(x::Group) = 1
iterate(x::Group) = (x, nothing)
iterate(x::Group, state) = nothing

function Group(individuals::Vector{T1}, policy::T2, pr::T3, meeting_days::Vector{Int}) where
    {T1<:Individual,T2<:Policy,T3<:Real}

    Group{T1,T2,T3}(string(uuid4()), individuals, policy, pr, meeting_days)
end

function now(group::Group)
    days = now.(group.individuals)
    if length(unique(days)) > 1
        throw(InexactError("Group is not time-synced"))
    end
    days[1]
end

test_and_isolate!(group::Group) = test_and_isolate!(group.policy, group)

function meet!(group::Group)
    n = length(group.individuals)
    if n <= 1
        return
    end
    if mod(now(group), 7) in group.meeting_days
        for i = 1:(n - 1), j = (i + 1):n
            a = group.individuals[i]
            b = group.individuals[j]
            if !is_isolating(a) & !is_isolating(b)
                if rand(Distributions.Bernoulli(group.meeting_probability))
                    meet!(a, b)
                end
            end
        end
    end
end

is_member(indv::Individual, g::Group) = indv in g.individuals

function get_adjacency_matrix(group::Group)
    n = length(group.individuals)
    A = Matrix{Float64}(undef, n, n)
    for i = 1:n
        A[i, i] = Inf64
    end
    for i = 1:(n - 1), j = (i + 1):n
        a = group.individuals[i]
        b = group.individuals[j]
        A[i, j] = 0.0
        if is_member(a, group) & is_member(b, group)
            A[i, j] += group.meeting_probability
        end
        A[j, i] = A[i, j]
    end
    return A
end
