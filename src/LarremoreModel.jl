struct LarremoreModel{T} <: DiseaseModel
    frac_symptomatic::T

    l10vl_onset::T
    day_onset_min::T
    day_onset_max::T

    l10vl_peak_min::T
    l10vl_peak_max::T
    peak_delay_max::T
	peak_delay_shape::T

	symptom_delay_min::T
    symptom_delay_max::T

    l10vl_clearance::T
	clearance_delay_min::T
    clearance_delay_max::T

	l10vl_infectivity_slope::T
end
function LarremoreModel(
    l10vl_infectivity_slope; # this is so completely arbitrary that a default makes no sense!
    frac_symptomatic = 0.75,
    l10vl_onset = 1e3, day_onset_min = 2.5, day_onset_max = 3.5,
    l10vl_peak_min = 7.0, l10vl_peak_max = 11.0, peak_delay_max = 3.0, peak_delay_shape = 1.5,
	symptom_delay_min = 0.0, symptom_delay_max = 3.0,
    l10vl_clearance = 6.0, clearance_delay_min = 4.0, clearance_delay_max = 9.0
)
    LarremoreModel{Float64}(
        frac_symptomatic,
        l10vl_onset, day_onset_min, day_onset_max,
        l10vl_peak_min, l10vl_peak_max, peak_delay_max, peak_delay_shape,
        symptom_delay_min, symptom_delay_max,
        l10vl_clearance, clearance_delay_min, clearance_delay_max,
        l10vl_infectivity_slope
    )
end

function sample(dm::LarremoreModel{T}) where {T<:Real}

    has_symptoms = rand(Distributions.Bernoulli(dm.frac_symptomatic))

    l10vl0 = 0
    l10vl_onset = dm.l10vl_onset
    l10vl_peak = rand(Distributions.Uniform(dm.l10vl_peak_min, dm.l10vl_peak_max))
    l10vl_clearance = dm.l10vl_clearance
    l10final = 0

    t0 = rand(Distributions.Uniform(-16/24, 8/24)) # sample infection time
    t1 = t0 + rand(Distributions.Uniform(dm.day_onset_min, dm.day_onset_max))
    tpeak = t1 + rand(Distributions.Uniform(dm.day_onset_min, dm.day_onset_max))
    tsymptoms = tpeak + (has_symptoms ? rand(Distributions.Uniform(dm.symptom_delay_min, dm.symptom_delay_max)) : 0)
    tclearance = tsymptoms + rand(Distributions.Uniform(dm.clearance_delay_min, dm.clearance_delay_max))
    clearance_slope = (l10vl_clearance - l10vl_peak) / (tclearance - tpeak)
    tfinal = tclearance - l10vl_clearance / clearance_slope

    t     = [t0,          t1,      tpeak,      tclearance, tfinal]
    l10vl = [ 1, l10vl_onset, l10vl_peak, l10vl_clearance,      0]
    approxfun = Interpolations.LinearInterpolation(t, l10vl; extrapolation_bc = Interpolations.Flat())
    tout = collect(0:ceil(tfinal)) .+ 7.5/24
    vlout = 10 .^ approxfun(tout)
    symptomatic = has_symptoms ? tout .>= tsymptoms : falses(length(tout))

    DiseaseTrajectory{T}(vlout, symptomatic)
end


function get_infection_probability(dm::LarremoreModel{T1}, dt::DiseaseTrajectory{T1}, day::T2) where
    {T1<:Real,T2<:Int}

    l10evl = max(0.0, log(10, get_viral_load(dt, day)) - dm.l10vl_clearance)
    min(1.0, max(0.0, dm.l10vl_infectivity_slope * l10evl))
end
