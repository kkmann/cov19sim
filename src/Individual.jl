mutable struct Individual{T1,T2,T3}
    uuid::UUID
    dm::T1
    dt::DiseaseTrajectory{T2}
    day_infected::T3
    symptom_probability::T2
    current_day::T3

    # logs etc
    isolation_log::Vector{Bool}
    isolation_stack::Vector{Bool}

    contacts_day::Vector{T3}
    contacts_id::Vector{UUID}
    contacts_infected::Vector{Bool}
    contacts_got_infected::Vector{Bool}

    test_log_day::Vector{T3}
    test_log_type::Vector{String}
    test_log_vl::Vector{T2}
    test_log_pr_pos::Vector{T2}
    test_log_result::Vector{Bool}
end

length(indv::Individual) = 1
iterate(indv::Individual) = (indv, nothing)
iterate(indv::Individual, state) = nothing

function Individual(dm::T1, symptom_probability::T2, isolation_timeframe::Int = 30) where {T1<:DiseaseModel,T2<:Real}
    Individual{T1,T2,Int}(
        uuid4(),
        dm, sample(dm), 2^30, symptom_probability, 0,
        falses(1), falses(isolation_timeframe),
        Vector{Int}(undef, 0), Vector{UUID}(undef, 0), Vector{Bool}(undef, 0), Vector{Bool}(undef, 0),
        Vector{Int}(undef, 0), Vector{String}(undef, 0), Vector{T2}(undef, 0), Vector{T2}(undef, 0), Vector{Bool}(undef, 0)
    )
end

function infect!(indv::Individual)
    if indv.day_infected == 2^30
        indv.day_infected = indv.current_day
    end
end

now(indv::Individual{T1,T2,T3}) where {T1<:DiseaseModel,T2<:Real,T3<:Int} = indv.current_day

get_viral_load(indv::Individual, day::T) where {T<:Int} = get_viral_load(indv.dt, day - indv.day_infected)
get_viral_load(indv::Individual) = get_viral_load(indv, indv.current_day)

is_symptomatic(indv::Individual, day::T) where {T<:Int} = is_symptomatic(indv.dt, day - indv.day_infected)
is_symptomatic(indv::Individual) = is_symptomatic(indv, indv.current_day)

function has_recovered(indv::Individual, day::T) where {T<:Int}
    has_recovered(indv.dt, day - indv.day_infected)
end
has_recovered(indv::Individual) = has_recovered(indv, indv.current_day)

get_infection_probability(indv::Individual, day::T) where {T<:Int} = get_infection_probability(indv.dm, indv.dt, day - indv.day_infected)
get_infection_probability(indv::Individual) = get_infection_probability(indv, indv.current_day)

is_infected(indv::Individual, day::T) where {T<:Int} = indv.day_infected <= day
is_infected(indv::Individual) = is_infected(indv, indv.current_day)

function is_isolating(indv::Individual, day::T) where {T<:Int}
    if day < 0
        return false
    end
    if day > indv.current_day
        throw(InexactError("day cannot lie in the future"))
    end
    return indv.isolation_log[day + 1]
end
is_isolating(indv::Individual) = is_isolating(indv, indv.current_day)

function isolate!(indv::Individual, isolation_pattern::Union{BitArray{1}, Vector{Bool}}; add_only = false)
    for i = 1:length(isolation_pattern)
        if add_only
            # only add new time
            indv.isolation_stack[i] = isolation_pattern[i] ? true : indv.isolation_stack[i]
        else
            # overwrite completely!
            indv.isolation_stack[i] = isolation_pattern[i] # remember how long to isolate for (-1 for including today!)
        end
    end
    # apply the first element immediately
    indv.isolation_log[indv.current_day + 1] = popat!(indv.isolation_stack, 1)
    # and fill up again
    push!(indv.isolation_stack, false)
end

isolate!(indv::Individual, duration::T) where {T<:Int} = isolate!(indv, trues(duration))

function isolate!(indv::Individual, from::T, to::T) where {T<:Int}

end

function log_contact!(a::Individual, b::Individual, a_infected_b::Bool, a_got_infected::Bool)
    push!(a.contacts_day, a.current_day)
    push!(a.contacts_id, b.uuid)
    push!(a.contacts_infected, a_infected_b)
    push!(a.contacts_got_infected, a_got_infected)
end

function log_test!(indv::Individual, type::String, result::Bool, pr_pos::T) where {T<:Real}
    push!(indv.test_log_day, indv.current_day)
    push!(indv.test_log_type, type)
    push!(indv.test_log_result, result)
    push!(indv.test_log_vl, get_viral_load(indv, indv.current_day))
    push!(indv.test_log_pr_pos, pr_pos)
end

function pcr_test_and_isolate!(indv::Individual, turnaround::Int, isolation::Int)
    # we assume perfect pcr testing
    pcr_positive = is_infected(indv)
    duration = turnaround + (pcr_positive ? isolation : 0)
    log_test!(indv, "pcr", pcr_positive, pcr_positive)
    isolate!(indv, duration)
    return pcr_positive
end

function step!(indv::Individual)
    # move the first element from the isolation stack (future behaviour) to
    # the isolation log (past behaviour)
    push!(indv.isolation_log, popat!(indv.isolation_stack, 1))
    push!(indv.isolation_stack, false) # add no isolation at the end of the stack
    # increase day counter by one
    indv.current_day += 1
end
function steps!(indv::Individual, n::Int)
    for i = 1:n
        step!(indv)
    end
end

function meet!(a::Individual, b::Individual)
    if a.current_day != b.current_day
        throw(InexactError("a and b must be time-synced"))
    end
    if (is_isolating(a) | is_isolating(b))
        throw(InexactError("neither a nor b can be isolating"))
    end
    if (is_infected(a) & is_infected(b)) | (!is_infected(a) & !is_infected(b))
        log_contact!(a, b, false, false)
        log_contact!(b, a, false, false)
        return # nothing to do, both infected
    end
    # at least one infected!
    if is_infected(a)
        infected = a
        noninfected = b
    else
        infected = b
        noninfected = a
    end
    if rand(Distributions.Bernoulli(get_infection_probability(infected)))
        infect!(noninfected)
        log_contact!(noninfected, infected, false, true)
        log_contact!(infected, noninfected, true, false)
    else
        log_contact!(noninfected, infected, false, false)
        log_contact!(infected, noninfected, false, false)
    end
end

function get_contact_log(indv::Individual)
    DataFrames.DataFrame(
        uuid_a = string(indv.uuid),
        uuid_b = string.(indv.contacts_id),
        day = indv.contacts_day,
        got_infected = indv.contacts_got_infected,
        infected_other = indv.contacts_infected
    )
end
get_contact_logs(indvs::Vector{T}) where {T<:Individual} = vcat(get_contact_log.(indvs)...)

function get_test_log(indv::Individual)
    DataFrames.DataFrame(
        uuid = string(indv.uuid),
        day = indv.test_log_day,
        type = indv.test_log_type,
        result = indv.test_log_result,
        viral_load = indv.test_log_vl,
        probability_positive = indv.test_log_pr_pos
    )
end
get_test_logs(indvs::Vector{T}) where {T<:Individual} = vcat(get_test_log.(indvs)...)

function get_status_log(indv::Individual)
    days = collect(0:(indv.current_day))
    function status(infected, infection_probability, symptomatic, recovered, isolated)
        if recovered
            return "recovered"
        elseif isolated
            return "isolated"
        elseif !infected
            return "susceptible"
        elseif (infection_probability > 0) & symptomatic
            return "infectious, symptomatic"
        elseif (infection_probability > 0) & !symptomatic
            return "infectious, asymptomatic"
        elseif infection_probability <= 0
            return "non-infectious"
        else
            return "NA"
        end
    end
    res = DataFrames.DataFrame(
        uuid = string(indv.uuid),
        day = days,
        infected = is_infected.(indv, days),
        viral_load = get_viral_load.(indv, days),
        symptomatic = is_symptomatic.(indv, days),
        infection_probability = get_infection_probability.(indv, days),
        recovered = has_recovered.(indv, days),
        isolated = is_isolating.(indv, days)
    )
    DataFrames.transform!(res,
        [:infected, :infection_probability, :symptomatic, :recovered, :isolated] =>
        ( (infected, infection_probability, symptomatic, recovered, isolated) ->
            status.(infected, infection_probability, symptomatic, recovered, isolated) ) =>
        :status
    )
    res
end
get_status_logs(indvs::Vector{T}) where {T<:Individual} = vcat(get_status_log.(indvs)...)
