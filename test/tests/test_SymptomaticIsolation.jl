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
    policy = SymptomaticIsolation(2, 8)
)
infect!.(school1.individuals[randperm(n_individuals(school1))[1:5]])
steps!(school1, 7*12)

school2 = sample_school(
    policy = SymptomaticIsolation(2, 19)
)
infect!.(school2.individuals[randperm(n_individuals(school2))[1:5]])
steps!(school2, 7*12)

println(evaluate([school1, school2]))
