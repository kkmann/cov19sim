abstract type Population end

length(x::Population) = 1
iterate(x::Population) = (x, nothing)
iterate(x::Population, state) = nothing

function step!(pop::Population)
    test_and_isolate!.(pop.groups)
    meet!.(pop.groups)
    step!.(pop.individuals)
end

function steps!(pop::Population, n::Int)
    for i = 1:n
        step!(pop)
    end
end

get_status_logs(pop::Population) = get_status_logs(pop.individuals)
get_contact_logs(pop::Population) = get_contact_logs(pop.individuals)
get_test_logs(pop::Population) = get_test_logs(pop.individuals)

function get_status_over_time(pop::Population)
    res = DataFramesMeta.@combine(
            DataFrames.groupby(get_status_logs(pop.individuals), [:day, :status]),
            count = length(:day)
        )
    res = DataFrames.unstack(res, :day, :status, :count)
    for col in names(res)
        res[!, col] .= Base.coalesce.(res[!, col], 0)
    end
    res
end

function get_contacts_over_time(pop::Population)
    DataFramesMeta.@combine(
        DataFrames.groupby(get_contact_logs(pop.individuals), :day),
        contacts = length(:uuid_a)
    )
end

function get_test_logs_over_time(pop::Population)
    DataFramesMeta.@combine(
        DataFrames.groupby(get_test_logs(pop.individuals), [:type, :day]),
        count = length(:uuid)
    )
end

function evaluate(pop::Population)
    n =  length(pop.individuals)
    n_days = pop.individuals[1].current_day
    logs = DataFramesMeta.@where(get_status_logs(pop), mod.(:day, 7) .<= 4)
    tmp = DataFramesMeta.@combine(
        DataFrames.groupby(
        DataFramesMeta.@where(
            DataFramesMeta.@combine(
                DataFrames.groupby(get_status_logs(pop.individuals), [:day, :status]),
                count = length(:day)
            ),
            (:status .== "infectious, asymptomatic") .| (:status .== "infectious, symptomatic")
        ), :day),
        count = sum(:count))
    infectious = Base.coalesce.(
        DataFrames.leftjoin(
            DataFrames.DataFrame(day = 0:n_days),
            tmp,
            on = :day
        )[!, :count], 0)
    contacts = Base.coalesce.(
        DataFrames.leftjoin(
            DataFrames.DataFrame(day = 0:n_days),
            get_contacts_over_time(pop),
            on = :day
        )[!, :contacts], 0) ./ n
    tests = DataFramesMeta.@where(get_test_logs_over_time(pop), :type == "pcr")
    DataFrames.DataFrame(
        id = string(uuid4()),
        proportion_infected = sum([is_infected(ind) for ind in pop.individuals]) / n,
        prop_workdays_missed = sum(logs[!, :isolated]) / DataFrames.nrow(logs),
        mean_daily_infectious = Statistics.mean(infectious),
        std_daily_infectious = Statistics.std(infectious),
        mean_daily_pp_contacts = Statistics.mean(contacts),
        mean_daily_pcr = Statistics.mean(tests[!, :count])
    )
end
