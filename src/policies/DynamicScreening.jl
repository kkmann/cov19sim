mutable struct DynamicScreening{T,PCR} <: Policy
    screening_test::T
    pcr_turnaround::Int
    isolation_duration::Int
    followup_duration::Int
    followup_weekdays::Vector{Int}
    fixed_isolation_weekdays::Vector{Int}
    screening_test_weekdays::Vector{Int}
    pcr_test::PCR
    buffer::Int
end



function DynamicScreening(screening_test::T;
    pcr_turnaround::Int = 2, isolation_duration::Int = 10,
    followup_duration::Int = 7, followup_weekdays::Vector{Int} = collect(0:4),
    fixed_isolation_weekdays::Vector{Int} = Int[],
    screening_test_weekdays::Vector{Int} = Int[],
    pcr_test::PCR = standard_pcr_test
) where {T<:Test,PCR<:Test}

    if type(pcr_test) != "pcr"
        throw(InexactError("pcr test must be called 'pcr'!"))
    end
    DynamicScreening{T,PCR}(screening_test, pcr_turnaround, isolation_duration, followup_duration,
        followup_weekdays, fixed_isolation_weekdays, screening_test_weekdays, pcr_test,
        0
    )
end



function test_and_isolate!(pol::DynamicScreening{T1}, gr::Group) where{T1<:Test}
    # apply flat isolation
    if mod(now(gr), 7) in pol.fixed_isolation_weekdays
        for x in gr.individuals
            isolate!(x, 1)
        end
    end
    ANY_PCR_POSITIVE = false
    # check for symptoms first
    ANY_SYMPTOMATIC = false
    for x in gr.individuals
        # check for new symptoms
        if is_newly_symptomatic(x)
            ANY_SYMPTOMATIC = true
            if is_positive(conduct_test!(pol.pcr_test, x))
                ANY_PCR_POSITIVE = true
                isolate!(x, pol.isolation_duration) # full duration
            else
                isolate!(x, pol.pcr_turnaround) # only PCR turnaround time
            end
        end
    end
    # handle screening tests
    # only pop one day from the follow-up buffer if it is the correct weekday
    # screen rest of group if anyone became symptomatic
    SCREEN_TODAY = ANY_SYMPTOMATIC
    if (mod(now(gr), 7) in pol.screening_test_weekdays)
        SCREEN_TODAY = true
    end
    if (pol.buffer > 0) & (mod(now(gr), 7) in pol.followup_weekdays)
        SCREEN_TODAY = true
    end
    ANY_SCREENING_POSITIVE = false
    if SCREEN_TODAY
        if pol.buffer > 0 # could also be due to regular screening!
            pol.buffer -= 1
        end
        for x in gr.individuals
            is_isolating(x) ? continue : nothing # skip already isolating individuals
            is_pcr_positive(x) ? continue : nothing # skip already pcr positive (already were isolated and are immune)
            if is_positive(conduct_test!(pol.screening_test, x))
                ANY_SCREENING_POSITIVE = true
                if is_positive(conduct_test!(pol.pcr_test, x)) # PCR follow-up
                    ANY_PCR_POSITIVE = true
                    isolate!(x, pol.isolation_duration) # full duration
                else
                    isolate!(x, pol.pcr_turnaround) # only PCR turnaround time
                end
            end
        end
    end
    # handle extension of dynamic screening policy: if anyone is newly symptomatic or tested positive,
    # extend dynamic testing period accordingly
    if ANY_SCREENING_POSITIVE | ANY_SYMPTOMATIC
        # -1 since we already screened once (today)
        buffer = (ANY_PCR_POSITIVE ? pol.followup_duration : pol.pcr_turnaround) - 1
    end
end
