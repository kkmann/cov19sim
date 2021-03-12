dm = LarremoreModel(
    0.75,
    3.0, 2.5, 3.5, # onset
    7.0, 11.0, 3.0, 1.5, # peak
    0.0, 3.0, # symptoms
    6.0, 4.0, 9.0, # clearance
    .25 # infectivtiy
)

indv1 = Individual(dm, 0.01; a = 2/30, b = 1/30)

println([Individual(dm, 0.01; a = 2/3, b = 1/3).compliance for i = 1:10])

@test indv1.current_day == 0

im = LogPropInfectionModel(.01, 1.0)
