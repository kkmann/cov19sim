seed!(42)

test = FixedTest("lfd", .5, .997)

function f(;
		gamma = 0.25,
		schooldays = collect(0:4),
		policy = DoNothing(),
		pr_symptoms = 0.01
	)

	dm = LarremoreModel(gamma; frac_symptomatic = 0.5)
	school = ThreeLevelPopulation(
		policy_bubble = policy,
		meeting_days = schooldays,
		disease_model = dm,
		pr_unrelated_symptoms = pr_symptoms
	)
	infect!.(school.individuals[randperm(n_individuals(school))[1:5]])
	steps!(school, 7*6)
    school
end

show(evaluate( [
	f(policy = DynamicScreening(test)),
	f(policy = DynamicScreening(test, screening_test_weekdays = [0, 2]))
] ))
