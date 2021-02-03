seed!(42)

function sample_school(;
    gamma = 0.05,
    schooldays = collect(0:4),
    policy = DoNothing(),
    pr_symptoms = 0.00
)
	dm = LarremoreModel(
		gamma;
	    frac_symptomatic = 0.5
	)
    ThreeLevelPopulation(
	    policy_bubble = policy,
	    meeting_days = schooldays,
	    disease_model = dm,
	    pr_unrelated_symptoms = pr_symptoms
    )
end

school1 = sample_school(
    policy = SymptomaticIsolation(0, true)
)
infect!.(school1.individuals[randperm(n_individuals(school1))[1:5]])
steps!(school1, 7*12)

school2 = sample_school(
    policy = SymptomaticIsolation(7, true)
)
infect!.(school2.individuals[randperm(n_individuals(school2))[1:5]])
steps!(school2, 7*12)

school3 = sample_school(
    policy = SymptomaticIsolation(21, true)
)
infect!.(school3.individuals[randperm(n_individuals(school3))[1:5]])
steps!(school3, 7*12)

println(evaluate([school1, school2, school3]))
