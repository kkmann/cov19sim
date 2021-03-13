dm = LarremoreModel(0.033)

individuals = [Individual(dm, 0.01) for i in 1:10]
infect!.(individuals)
lfd = LogRegTest("lfd", 0.76, -3.81, .998)
for i in 1:30
    step!.(individuals)
    conduct_test!.(lfd, individuals)
end

tmp = get_test_logs(individuals)



# check that autoregressive components works!
individual = Individual(dm, 0.0)
infect!(individual)
# fix the viral load
individual.dt.vl .= 10.0
lfd = FixedTest("lfd", .5, 1.0; ar_window = 2, ar_coefficient = 0.75)
@test get_probability_positive(lfd, individual) == 0.5
step!(individual)
@test get_probability_positive(lfd, individual) == 0.5
conduct_test!(lfd, individual)
individual.test_log_result[1] = true # guarantee that the test is positive
@test get_probability_positive(lfd, individual) == 0.75 + 0.25*0.5
step!(individual)
step!(individual)
@test get_probability_positive(lfd, individual) == 0.75 + 0.25*0.5
step!(individual)
@test get_probability_positive(lfd, individual) == 0.5



pcr = FixedTest("pcr", 0.975, 1.0; lod = 300.0)
@test get_probability_positive(pcr, 150.0) == 0.0
@test get_probability_positive(pcr, 300.0) == 0.975
@test get_probability_positive(pcr, 1e11) == 0.975


# check that random effects for logreg test work
lfd = LogRegTest("lfd", 0.76, -3.81, .998; ranef = 0.0)
@test get_probability_positive(lfd, 1e7; u = 1.0) ==
    get_probability_positive(lfd, 1e7; u = 3.0)

lfd = LogRegTest("lfd", 0.76, -3.81, .998; ranef = 1.0)
@test get_probability_positive(lfd, 1e7; u = 1.0) <
    get_probability_positive(lfd, 1e7; u = 3.0)
@test Individual(dm, 0.0).u_sensitivity != 0.0
