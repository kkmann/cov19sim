school = ThreeLevelPopulation(
    policy_bubble = SymptomaticIsolation(14, true)
)
infect!.(school.individuals[1:5])
steps!(school, 80)
#@test all([ind.current_day for ind in school.individuals[1]] .== 80)
#tbl_status = get_status_logs(school)
#show(get_status_over_time(school))
#show(evaluate(school))

import Distributed.pmap, Distributed.addprocs, Distributed.@everywhere
addprocs(12)
@everywhere begin

using cov19sim
function f(i)
    school = ThreeLevelPopulation(policy_bubble = SymptomaticIsolation(14, true))
    infect!.(school.individuals[1])
    steps!(school, 80)
    evaluate(school)
end

end

@time show(vcat(pmap(f, 1:500)...))

tmp = @time vcat(pmap(f, 1:1000)...)
tmp = @time vcat(pmap(f, 1:2000)...)


tmp = @time vcat(map(f, 1:1000)...)
