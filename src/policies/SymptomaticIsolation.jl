struct SymptomaticIsolation{T,PCR} <: Policy
    pcr_turnaround::Int
    isolation_duration::Int
    screening_test_weekdays::Vector{Int}
    fixed_isolation_weekdays::Vector{Int}
    screening_test::T
    pcr_test::PCR
end
function SymptomaticIsolation(
        screening_test::T;
        pcr_turnaround::Int = 2,
        isolation_duration::Int = 10,
        screening_test_weekdays::Vector{Int} = Int[],
        fixed_isolation_weekdays::Vector{Int} = Int[],
        pcr_test::PCR = standard_pcr_test
    ) where {T<:Test,PCR<:Test}

    if type(pcr_test) != "pcr"
        throw(InexactError("pcr test must be called 'pcr'!"))
    end
    SymptomaticIsolation{T,PCR}(pcr_turnaround, isolation_duration, screening_test_weekdays, fixed_isolation_weekdays, screening_test, pcr_test)
end



function test_and_isolate!(pol::SymptomaticIsolation, gr::Group)
    # apply fixed-isolation
    if mod(now(gr), 7) in pol.fixed_isolation_weekdays
        for x in gr.individuals
            isolate!(x, 1)
        end
    end
    any_symptomatic = false
    any_screening_test_positive = false
    any_pcr_positive = false
    for x in gr.individuals
        # check symptoms first
        if is_newly_symptomatic(x) # check for new symptoms
            any_symptomatic = true
            pcr_test_positive = is_positive(conduct_test!(pol.pcr_test, x))
            if pcr_test_positive
                any_pcr_positive = true
                isolate!(x, pol.isolation_duration) # full duration
            else
                isolate!(x, pol.pcr_turnaround) # only PCR turnaround time
            end
        end
        # if screening weekday, apply screening tests to all non-isolating (excl. symtomatic)
        # skip individuals who were pcr positive previously (already underwent isolation)
        if mod(now(gr), 7) in pol.screening_test_weekdays
            if !is_isolating(x) & !is_newly_symptomatic(x) & !is_pcr_positive(x)
                screening_test_positive = is_positive(conduct_test!(pol.screening_test, x))
                if screening_test_positive # PCR-test + immediately isolate individuals (if positive)
                    any_screening_test_positive = true
                    pcr_test_positive = is_positive(conduct_test!(pol.pcr_test, x))
                    if pcr_test_positive
                        any_pcr_positive = true
                        isolate!(x, pol.isolation_duration) # full duration
                    else
                        isolate!(x, pol.pcr_turnaround) # only PCR turnaround time
                    end
                end
            end
        end
    end
    # handle bubble isolation
    if any_screening_test_positive
        if any_pcr_positive
            # all go into isolaion for full duration immediately
            for x in gr.individuals
                isolate!(x, pol.isolation_duration) # overwrites existing isolation
            end
        else
            # all go into isolaion for turnaroudn time only
            for x in gr.individuals
                isolate!(x, pol.pcr_turnaround) # overwrites existing isolation
            end
        end
    else
        # no screening positives, only go into isolation when positive PCR comes back
        if any_pcr_positive
            new_pattern = trues(pol.isolation_duration)
            new_pattern[1:(pol.pcr_turnaround)] .= false
            for x in gr.individuals
                isolate!(x, new_pattern; add_only = true) # do no overwrite non-isolation within turnaround!
            end
        end
    end
end
