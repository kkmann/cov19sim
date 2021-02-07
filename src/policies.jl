# generic methods
function test_and_isolate!(pol::T1, gr::Group) where {T1<:Policy}
    throw(MethodError("not implemented"))
end

length(x::Policy) = 1
iterate(x::Policy) = (x, nothing)
iterate(x::Policy, state) = nothing



struct DoNothing <: Policy end
test_and_isolate!(pol::DoNothing, gr::Group) = return



struct SymptomaticIsolation <: Policy
    pcr_turnaround::Int
    isolation_duration::Int
end

function test_and_isolate!(pol::SymptomaticIsolation, gr::Group)
    # screen all individuals, isolate all symptomatic and conduct pcr follow-up
    any_new_pcr_positive = false
    for x in gr.individuals
        if is_symptomatic(x)
            pcr_positive = pcr_test_and_isolate!(x, pol.pcr_turnaround,  pol.isolation_duration)
            any_new_pcr_positive = pcr_positive ? true : any_new_pcr_positive
        end
    end
    println(now(gr))
    println(any_new_pcr_positive)
    # adjust the isolation time for everyone
    if any_new_pcr_positive
        new_pattern = trues(pol.pcr_turnaround + pol.isolation_duration)
        new_pattern[1:(pol.pcr_turnaround)] .= false
        for x in gr.individuals
            isolate!(x, new_pattern; add_only = true) # we do not want to 'unisolate' if already isolating for other reasons
        end
    end
end



struct RegularScreening{T} <: Policy
    screening_test::T
    pcr::Bool
    pcr_turnaround::Int
    isolation_duration::Int
    isolate_all::Bool
    test_weekdays::Vector{Int}
end

function RegularScreening(screening_test::T;
    pcr::Bool = true, pcr_turnaround::Int = 2, isolation_duration::Int = 8,
    isolate_all::Bool = true, test_weekdays::Vector{Int} = collect(0:4)
) where {T<:Test}
    RegularScreening{T}(screening_test, pcr, pcr_turnaround, isolation_duration, isolate_all, test_weekdays)
end

function isolate!(pol::RegularScreening{T}, indv::Individual) where {T<:Test}
    if pol.pcr
        pcr_test_and_isolate!(indv, pol.pcr_turnaround, pol.isolation_duration)
    else
        isolate!(indv, pol.isolation_duration)
    end
end

function test_and_isolate!(pol::RegularScreening{T1}, gr::Group) where{T1<:Test}
    to_isolate = Vector{typeof(gr.individuals[1])}()
    triggered = false
    for x in gr.individuals
        if is_isolating(x)
            continue # check next
        end
        isolate = is_symptomatic(x) | triggered
        if !isolate & (mod(now(x), 7) in pol.test_weekdays) # screening test?
            isolate = conduct_test!(pol.screening_test, x)
        end
        if !isolate
            push!(to_isolate, x) # add to list of individuals who need to isolate if something comes up
        else
            isolate!(pol, x)
            triggered = true
        end
    end
    if triggered
        for x in to_isolate
            isolate!(pol, x)
        end
    end
end



# stateful! each group needs its own instance
mutable struct DynamicScreening{T} <: Policy
    screening_test::T
    pcr_turnaround::Int
    isolation_duration::Int
    followup_duration::Int
    followup_weekdays::Vector{Int}
    buffer::Int
    buffer_t::Int
end

function DynamicScreening(screening_test::T;
    pcr_turnaround::Int = 2, isolation_duration::Int = 8,
    followup_duration::Int = 7, followup_weekdays::Vector{Int} = collect(0:4)
) where {T<:Test}
    DynamicScreening{T}(screening_test, pcr_turnaround, isolation_duration, followup_duration,
        followup_weekdays, 0, 0
    )
end

function test_and_isolate!(pol::DynamicScreening{T1}, gr::Group) where{T1<:Test}
    screen_today = (pol.buffer > 0) & (mod(now(gr), 7) in pol.followup_weekdays)
    # update buffer, make sure we are in sync with group timeline
    if pol.buffer_t <= now(gr)
        pol.buffer_t += 1
        # check that the weekday is a followup day!
        if screen_today
            pol.buffer -= 1
        end
    end
    to_screen = Vector{typeof(gr.individuals[1])}() # buffer for individuals to be screened that day
    followup_extension = 0 # days to put on follow up
    for x in gr.individuals
        is_isolating(x) ? continue : nothing # skip isolating individuals
        # detect based on symptoms and screening test
        positive = is_symptomatic(x)
        if !positive & screen_today
            # check if screening tests comes back positive
            positive = conduct_test!(pol.screening_test, x)
        end
        # follow up with pcr and dermine extension of the dynamic testing scheme
        if positive
             pcr_positive = pcr_test_and_isolate!(x, pol.pcr_turnaround, pol.isolation_duration)
             # compute extension of the follow up regime based on this individual; -1 since testing starts today (incl)
             extension = max(0, (pcr_positive ? pol.followup_duration : pol.pcr_turnaround) - 1)
             # see if it is longer than already needed, update
             followup_extension = max(followup_extension, extension)
        else
            # negative, add to list of individuals who need to be screened if something pops up later
            push!(to_screen, x)
        end
    end
    # screen remaining individuals if positive extension (and correct weekday)
    if (followup_extension > 0) & screen_today
        for x in to_screen
            screening_positive = conduct_test!(pol.screening_test, x)
            if screening_positive
                # handle follow up if screening turns out positive
                pcr_positive =  pcr_test_and_isolate!(x, pol.pcr_turnaround, pol.isolation_duration)
                extension = max(0, (pcr_positive ? pol.followup_duration : pol.pcr_turnaround) - 1)
                followup_extension = max(followup_extension, extension)
            end
        end
    end
    # update dynamic screening buffer
    pol.buffer = max(pol.buffer, followup_extension)
end
