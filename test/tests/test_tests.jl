dm = LarremoreModel(0.033)

individuals = [Individual(dm, 0.01) for i in 1:10]
infect!.(individuals)
lfd2 = LogRegTest("lfd-2", 0.76, .5, -3.81, .997)
for i in 1:30
    step!.(individuals)
    conduct_test!.(lfd2, individuals)
end

tmp = get_test_logs(individuals)
