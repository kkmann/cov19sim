dm = LarremoreModel(0.033)

individuals = [Individual(dm, 0.01) for i in 1:10]
infect!.(individuals)
lfd1 = LogPropTest("lfd-1", 5.0, 1/12, .997)
lfd2 = LogGLMTest("lfd-2", 1e5, 2.0, .5, -15.0, .997)
for i in 1:30
    step!.(individuals)
    conduct_test!.(lfd1, individuals)
    conduct_test!.(lfd2, individuals)
end

show( get_test_logs(individuals) )
