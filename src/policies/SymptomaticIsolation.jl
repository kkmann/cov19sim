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
    # adjust the isolation time for everyone
    if any_new_pcr_positive
        new_pattern = trues(pol.pcr_turnaround + pol.isolation_duration)
        new_pattern[1:(pol.pcr_turnaround)] .= false
        for x in gr.individuals
            isolate!(x, new_pattern; add_only = true) # we do not want to 'unisolate' if already isolating for other reasons
        end
    end
end
