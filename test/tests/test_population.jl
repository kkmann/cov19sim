seed!(42)

dm = LarremoreModel(.05)
school = ThreeLevelPopulation(
    policy_bubble = SymptomaticIsolation(
        FixedTest("lfd", .5, .997);
        screening_test_weekdays = [0]
    ),
    meeting_days = collect(0:4),
    disease_model = dm;
    a = .66, b = .34
)
idx = randperm(n_individuals(school))[1:5]
infect!.(school.individuals[idx])
steps!(school, 7*12)

A = get_adjacency_matrix(school)
@test minimum(A) >= 0

resample(school)

println(evaluate(school; tests = ["pcr", "lfd"]))
