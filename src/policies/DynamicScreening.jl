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
