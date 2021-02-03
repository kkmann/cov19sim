seed!(42)

lfd_test = LogPropTest(
	"lfd",
	5.0, # log-10 VL LOD
	1/12, # log-10 excess VL proportionality constant for sensitivity
	.997 # specificity
)

function eval(;
    gamma = 0.05,
    schooldays = collect(0:4),
    policy = DoNothing(),
	pr_symptoms = 0.00
)
	dm = LarremoreModel(gamma; frac_symptomatic = 0.5)
	school = ThreeLevelPopulation(
		policy_bubble = policy,
		meeting_days = schooldays,
		disease_model = dm,
		pr_unrelated_symptoms = pr_symptoms
	)
	infect!.(school.individuals[randperm(n_individuals(school))[1:5]])
	steps!(school, 7*12)
	return evaluate(school)
end

show(vcat(
	eval(policy = RegularScreening(lfd_test; test_weekdays = collect(0:4))),
	eval(policy = RegularScreening(lfd_test; pcr = false, test_weekdays = collect(0:4))),
	eval(policy = RegularScreening(lfd_test; test_weekdays = [0])),
	eval(policy = RegularScreening(lfd_test; test_weekdays = [0]), gamma = 0.1)
))
