struct SplitGroup{T,PCR} <: Policy
    a_days::Vector{Int}
    b_days::Vector{Int}
    pcr_turnaround::Int
    isolation_duration::Int
    screening_test::T
    screening_test_weekdays::Vector{Int}
    pcr_test::PCR
end



function SplitGroup(
    screening_test::T;
    a_days = [0, 1],
    b_days = [2, 3],
    pcr_turnaround::Int = 2,
    isolation_duration::Int = 10,
    screening_test_weekdays::Vector{Int} = Int[],
    pcr_test::PCR = standard_pcr_test
) where {T<:Test,PCR<:Test}

    if type(pcr_test) != "pcr"
        throw(InexactError("pcr test must be called 'pcr'!"))
    end
    SplitGroup{T,PCR}(a_days, b_days, pcr_turnaround, isolation_duration, screening_test, screening_test_weekdays, pcr_test)
end

is_a(i::Int, pol::SplitGroup) = mod(i, 2) == 1

function test_and_isolate!(pol::SplitGroup, gr::Group)
    n = length(gr.individuals)
    today = now(gr)
    # isolate everyone according to alternating scheme
    for i = 1:n
        indv = gr.individuals[i]
        if is_a(i, pol) & !(mod(today, 7) in pol.a_days) & !is_isolating(indv)
            isolate!(indv, 1)
        end
        if !is_a(i, pol) & !(mod(today, 7) in pol.b_days) & !is_isolating(indv)
            isolate!(indv, 1)
        end
    end
    # check for new symptoms / screeing tests in sub groups
    ANY_SYMPTOMATIC_A = false
    ANY_SCREENING_POSITIIVE_A = false
    ANY_PCR_POSITIIVE_A = false
    ANY_SYMPTOMATIC_B = false
    ANY_SCREENING_POSITIIVE_B = false
    ANY_PCR_POSITIIVE_B = false
    for i in 1:n
        indv = gr.individuals[i]
        if is_symptomatic(indv)
            if is_a(i, pol)
                ANY_SYMPTOMATIC_A = true
            else
                ANY_SYMPTOMATIC_B = true
            end
            pcr_positive = conduct_test!(pol.pcr_test, indv)
            if pcr_positive
                if is_a(i, pol)
                    ANY_PCR_POSITIIVE_A = true
                else
                    ANY_PCR_POSITIIVE_B = true
                end
                isolate!(indv, pol.isolation_duration)
            else
                isolate!(indv, pol.pcr_turnaround)
            end
        else
            if mod(now(gr), 7) in pol.screening_test_weekdays
                if conduct_test!(pol.screening_test, indv)
                    if is_a(i, pol)
                        ANY_SCREENING_POSITIIVE_A = true
                    else
                        ANY_SCREENING_POSITIIVE_B = true
                    end
                    pcr_positive = conduct_test!(pol.pcr_test, indv)
                    if pcr_positive
                        if is_a(i, pol)
                            ANY_PCR_POSITIIVE_A = true
                        else
                            ANY_PCR_POSITIIVE_B = true
                        end
                        isolate!(indv, pol.isolation_duration)
                    else
                        isolate!(indv, pol.pcr_turnaround)
                    end
                end
            end
        end
    end
    # handle subgroup isolation
    for i in 1:n
        indv = gr.individuals[i]
        if is_a(i, pol)
            if ANY_SCREENING_POSITIIVE_A
                if ANY_PCR_POSITIIVE_A
                    isolate!(indv, pol.isolation_duration)
                else
                    isolate!(indv, pol.pcr_turnaround)
                end
            else
                if ANY_PCR_POSITIIVE_A
                    new_pattern = trues(pol.isolation_duration)
                    new_pattern[1:(pol.pcr_turnaround)] .= false
                    isolate!(indv, new_pattern; add_only = true)
                end
            end
        else # b
            if ANY_SCREENING_POSITIIVE_B
                if ANY_PCR_POSITIIVE_B
                    isolate!(indv, pol.isolation_duration)
                else
                    isolate!(indv, pol.pcr_turnaround)
                end
            else
                if ANY_PCR_POSITIIVE_B
                    new_pattern = trues(pol.isolation_duration)
                    new_pattern[1:(pol.pcr_turnaround)] .= false
                    isolate!(indv, new_pattern; add_only = true)
                end
            end
        end
    end
end
