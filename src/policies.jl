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
    isolation_duration::Int
    isolate_all::Bool
end

function test_and_isolate!(pol::SymptomaticIsolation, gr::Group)
    for x in gr.individuals
        is_isolating(x) ? continue : nothing # skip isolating individuals
        if is_symptomatic(x)
            if pol.isolate_all
                for xx in gr.individuals
                    !is_isolating(xx) ? isolate!(xx, pol.isolation_duration) : nothing
                end
            else
                isolate!(x, pol.isolation_duration)
            end
            break # no need to check the others!
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
    pcr::Bool = true, pcr_turnaround::Int = 2, isolation_duration::Int = 14, isolate_all::Bool = true, test_weekdays::Vector{Int} = collect(0:4)
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
    for x in gr.individuals
        if is_isolating(x)
            continue # check next
        end
        isolate = is_symptomatic(x)
        if !isolate & (now(x) in pol.test_weekdays) # screening test?
            global isolate = conduct_test!(pol.screening_test, x)
        end
        if !isolate
            continue # nothing to do
        end
        if pol.isolate_all
            for xx in gr.individuals
                isolate!(gr, x)
            end
        else
            isolate!(gr, x)
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
    pcr_turnaround::Int = 2, isolation_duration::Int = 14, followup_duration::Int = 7, followup_weekdays::Vector{Int} = collect(0:4)
) where {T<:Test}
    DynamicScreening{T}(screening_test, pcr_turnaround, isolation_duration, followup_duration,
        followup_weekdays, 0, 0
    )
end

function test_and_isolate!(pol::DynamicScreening{T1}, gr::Group) where{T1<:Test}
    # update buffer, make sure we are in sync with group timeline
    if pol.buffer_t <= now(gr.individuals[1])
        pol.buffer_t += 1
        # check that the weekday is a followup day!
        if (pol.buffer > 0) & (now(gr.individuals[1]) in pol.followup_weekdays)
            pol.buffer -= 1
        end
    end
    to_screen = Vector{typeof(gr.individuals[1])}() # buffer for individuals to be screened that day
    followup_extension = 0 # days to put on follow up
    for x in gr.individuals
        is_isolating(x) ? continue : nothing # skip isolating individuals
        # detect based on symptoms and screening test
        positive = is_symptomatic(x)
        if !positive & (pol.buffer > 0) & (now(x) in pol.followup_weekdays)
            # check if screening tests comes back positive
            positive = conduct_test!(pol.screening_test, x)
        end
        # follow up with pcr and dermine extension of the dynamic testing scheme
        if positive
             pcr_positive = pcr_test_and_isolate!(x, pol.pcr_turnaround, pol.isolation_duration)
             # compute extension of the follow up regime based on this individual; -1 since testing starts today (incl)
             extension = max(0, (pcr_positive ? pol.followup_duration : pol.pcr_turnaround) - 1)
             # see if it is longer than already needed, update
             global followup_extension = max(followup_extension, extension)
        else
            # negative, add to list of individuals who need to be screened if something pops up later
            push!(to_screen, x)
        end
    end
    # screen remaining individuals if positive extension (and correct weekday)
    if (followup_extension > 0) & (now(gr.individuals[1]) in pol.followup_weekdays)
        for x in to_screen
            screening_positive = conduct_test!(pol.screening_test, x)
            if screening_positive
                # handle follow up if screening turns out positive
                pcr_positive =  pcr_test_and_isolate!(x, pol.pcr_turnaround, pol.isolation_duration)
                extension = max(0, (pcr_positive ? pol.followup_duration : pol.pcr_turnaround) - 1)
                global followup_extension = max(followup_extension, extension)
            end
        end
    end
    # add to dynamic screening buffer
    pol.buffer += followup_extension
end
