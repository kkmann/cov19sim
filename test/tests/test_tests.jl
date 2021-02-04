dm = LarremoreModel(0.033)

individuals = [Individual(dm, 0.01) for i in 1:10]
infect!.(individuals)
lfd1 = LogPropTest("lfd-1", 5.0, 1/12, .997)
lfd21 = LogGLMTest("lfd-2", 5.0, 1/12, .997)
for i in 1:30
    step!.(individuals)
    conduct_test!.(lfd1, individuals)
end

tmp = get_test_logs(individuals)
