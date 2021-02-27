mutable struct ThreeLevelPopulation{T1} <: Population
    individuals::Vector{T1}
    groups::Vector{Any}
    initial_parameters::Dict
end

function ThreeLevelPopulation(;
    n_per_bubble = 7,
    bubbles_per_class = 3,
    m_classes = 12,
    pr_meet_bubble = 0.9,
    pr_meet_class = 0.33,
    pr_meet_school = 1/251,
    policy_bubble = DoNothing(),
    policy_class = DoNothing(),
    policy_school = DoNothing(),
    meeting_days = collect(0:4),
    disease_model = LarremoreModel(0.033),
    pr_unrelated_symptoms = 0.01,
    a = Inf,
    b = 1.0
)
    initial_parameters = Dict(
        :n_per_bubble => n_per_bubble,
        :bubbles_per_class => bubbles_per_class,
        :m_classes => m_classes,
        :pr_meet_bubble => pr_meet_bubble,
        :pr_meet_class => pr_meet_class,
        :pr_meet_school => pr_meet_school,
        :policy_bubble => policy_bubble,
        :policy_class => policy_class,
        :policy_school => policy_school,
        :meeting_days => meeting_days,
        :disease_model => disease_model,
        :pr_unrelated_symptoms => pr_unrelated_symptoms,
        :a => a,
        :b => b
    )
    bubbles = []
    for i in 1:(bubbles_per_class*m_classes)
        bubble_individuals = [Individual(disease_model, pr_unrelated_symptoms; a = a, b = b) for j = 1:n_per_bubble]
        push!(bubbles, Group(bubble_individuals, deepcopy(policy_bubble), pr_meet_bubble, meeting_days))
    end
    classes = []
    for i in 1:m_classes
        class_individuals = vcat([bubbles[j].individuals for j = 1:length(bubbles) if (mod(j, m_classes) + 1 == i)]...)
        push!(classes, Group(class_individuals, deepcopy(policy_class), pr_meet_class, meeting_days))
    end
    individuals = vcat([classes[i].individuals for i = 1:length(classes)]...)
    school = Group(individuals, deepcopy(policy_school), pr_meet_school, meeting_days)
    ThreeLevelPopulation{typeof(individuals[1])}(individuals, vcat(bubbles, classes, school), initial_parameters)
end

function resample(population::ThreeLevelPopulation{I}) where {I<:Individual}
    ThreeLevelPopulation(; population.initial_parameters...)
end
