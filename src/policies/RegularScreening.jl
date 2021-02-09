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
