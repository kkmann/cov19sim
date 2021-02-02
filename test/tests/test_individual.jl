dm = LarremoreModel(
    0.75,
    3.0, 2.5, 3.5, # onset
    7.0, 11.0, 3.0, 1.5, # peak
    0.0, 3.0, # symptoms
    6.0, 4.0, 9.0, # clearance
    .25 # infectivtiy
)

indv1 = Individual(dm, 0.01)

@test indv1.current_day == 0
indv1.current_day = 1
println(indv1.current_day)
println(indv1.day_infected)
infect!(indv1)
println(indv1.day_infected)


indv1 = Individual(dm, 0.01)
indv2 = Individual(dm, 0.01)

println(get_contact_log(indv1))

infect!(indv1)
steps!(indv1, 5)
steps!(indv2, 5)
println(get_infection_probability(indv1, 5))
for i = 1:25
    meet!(indv1, indv2)
end

println(get_contact_log(indv1))


println(get_status_log(indv1))

log_test!(indv1, "bla", true)
println(get_test_log(indv1))


println(get_status_logs([indv1, indv2]))
