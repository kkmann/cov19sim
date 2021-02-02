dm = LarremoreModel(
    0.75,
    3.0, 2.5, 3.5, # onset
    7.0, 11.0, 3.0, 1.5, # peak
    0.0, 3.0, # symptoms
    6.0, 4.0, 9.0, # clearance
    .25 # infectivtiy
)

gr = Group(
    [Individual(dm, 0.01) for i = 1:10],
    SymptomaticIsolation(14, true),
    .33,
    collect(0:4)
)

meet!(gr)

println(get_contact_logs(gr.individuals))

steps!.(gr.individuals, 5)
meet!(gr)
println(get_contact_logs(gr.individuals))

test_and_isolate!(gr)
