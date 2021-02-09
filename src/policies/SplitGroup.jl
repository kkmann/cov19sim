struct SplitGroup <: Policy
    a_days::Vector{Int}
    b_days::Vector{Int}
    pcr_turnaround::Int
    isolation_duration::Int
end



function SplitGroup(;
    a_days = [0, 2],
    b_days = [1, 3],
    pcr_turnaround::Int = 2,
    isolation_duration::Int = 8
)
    SplitGroup(a_days, b_days, pcr_turnaround, isolation_duration)
end

is_a(i::Int, pol::SplitGroup) = mod(i, 2) == 1
is_b(i::Int, pol::SplitGroup) = !is_a(i, pol)

function test_and_isolate!(pol::SplitGroup, gr::Group)
    n = length(gr.individuals)
    today = now(gr)
    # isolate everyone according to alternating scheme
    for i = 1:n
        indv = gr.individuals[i]
        if is_a(i, pol) & !(mod(today, 7) in pol.a_days) & !is_isolating(indv)
            isolate!(indv, 1)
        end
        if is_b(i, pol) & !(mod(today, 7) in pol.b_days) & !is_isolating(indv)
            isolate!(indv, 1)
        end
    end
    # check for new symptoms in sub groups
    triggered_a = false
    to_isolate_a = Vector{typeof(gr.individuals[1])}()
    triggered_b = false
    to_isolate_b = Vector{typeof(gr.individuals[1])}()
    for i in 1:n
        indv = gr.individuals[i]
        if is_symptomatic(indv)
            if pcr_test_and_isolate!(indv, pol.pcr_turnaround, pol.isolation_duration)
                is_a(i, pol) ? triggered_a = true : nothing
                is_b(i, pol) ? triggered_b = true : nothing
            end
        end
    end
    if triggered_a | triggered_b
        # isoalte all others from day of positive result (2 days from now)
        new_pattern = trues(pol.pcr_turnaround + pol.isolation_duration)
        new_pattern[1:(pol.pcr_turnaround)] .= false
        for i in 1:n
            indv = gr.individuals[i]
            if (triggered_a & is_a(i, pol)) | (triggered_b & is_b(i, pol))
                isolate!(indv, new_pattern; add_only = true) # we do not want to 'unisolate' if already isolating for other reasons
            end
        end
    end
end
