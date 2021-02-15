seed!(42)

function sample_school(;
    gamma = 0.25,
    schooldays = collect(0:4),
    policy = DoNothing(),
    pr_symptoms = 0.00
)
	dm = LarremoreModel(gamma; frac_symptomatic = 0.5)
    ThreeLevelPopulation(
	    policy_bubble = policy,
	    meeting_days = schooldays,
	    disease_model = dm,
	    pr_unrelated_symptoms = pr_symptoms
    )
end

test = FixedTest("lfd", .5, .997)

school1 = sample_school(
    policy = SymptomaticIsolation(test)
)
infect!.(school1.individuals[randperm(n_individuals(school1))[1:5]])
steps!(school1, 7*6)

school2 = sample_school(
    policy = SymptomaticIsolation(test, fixed_isolation_weekdays = [3, 4])
)
infect!.(school2.individuals[randperm(n_individuals(school2))[1:5]])
steps!(school2, 7*6)

school3 = sample_school(
    policy = SymptomaticIsolation(test, screening_test_weekdays = [0, 3])
)
infect!.(school3.individuals[randperm(n_individuals(school3))[1:5]])
steps!(school3, 7*6)

school4 = sample_school(
    policy = SymptomaticIsolation(test);
    pr_symptoms = 0.05
)
infect!.(school4.individuals[randperm(n_individuals(school4))[1:5]])
steps!(school4, 7*6)

println(evaluate([school1, school2, school3, school4]))
