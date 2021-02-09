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

school = sample_school(
    policy = SplitGroup()
)
infect!.(school.individuals[randperm(n_individuals(school1))[1:5]])
steps!(school, 7*12)

show(evaluate(school))
