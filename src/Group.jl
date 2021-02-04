
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
    if mod(group.individuals[1].current_day, 7) in group.meeting_days
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
