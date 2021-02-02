dm = LarremoreModel(
    0.75,
    3.0, 2.5, 3.5, # onset
    7.0, 11.0, 3.0, 1.5, # peak
    0.0, 3.0, # symptoms
    6.0, 4.0, 9.0, # clearance
    .25 # infectivtiy
)

indv = Individual(dm, 0.01)
infect!(indv)
steps!(indv, 5)

LFD = LogPropTest("LFD", 5.0, 1/12, .997)

println(get_viral_load.(indv, 0:5))

for i = 1:10
    println(conduct_test(LFD, indv))
end

println(get_test_log(indv))
