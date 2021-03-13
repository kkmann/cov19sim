abstract type Population end

length(x::Population) = 1
iterate(x::Population) = (x, nothing)
iterate(x::Population, state) = nothing

n_individuals(pop::Population) = length(pop.individuals)

resample(population::Population) = throw(MethoError("not implemented"))

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

function get_mean_contacts_over_time(pop::Population)
    res = DataFrames.leftjoin(
        DataFrames.DataFrame(
            day = 0:(now(pop.individuals[1]) - 1)
            ),
        DataFramesMeta.@combine(
            DataFrames.groupby(get_contact_logs(pop.individuals), :day),
            contacts = length(:uuid_a)
        ),
        on = :day
    )
    # replace missing with 0
    res[!, :contacts] .= Base.coalesce.(res[!, :contacts], 0) ./ n_individuals(pop)
    res
end

function get_test_logs_over_time(pop::Population)
    tmp = DataFramesMeta.@combine(
        DataFrames.groupby(get_test_logs(pop.individuals), [:type, :day]),
        count = length(:uuid)
    )
    unique_types = unique(tmp[!,:type])
    if length(unique_types) == 0
        return DataFrames.DataFrame( # no tests
            day = Vector{Int}(),
            type = Vector{String}(),
            count = Vector{Int}()
        )
    else
        tbl_grid = DataFrames.DataFrame(
            gridexpand(
                0:(now(pop.individuals[1]) - 1),
                unique_types
            ),
            ["day", "type"]
        )
        res = DataFrames.leftjoin(tbl_grid, tmp, on = [:day, :type])
        res[!, :day] .= Int.(res[!, :day])
        res[!, :type] .= string.(res[!, :type])
        res[!, :count] .= Base.coalesce.(res[!, :count], 0)
        return res
    end
end

n_infected(pop::Population) = sum(is_infected.(pop.individuals))

function n_workdays_missed(pop::Population; workdays::Vector{Int} = collect(0:4))
    DataFrames.nrow(
        DataFramesMeta.@where(get_status_logs(pop), is_workday.(:day; workdays = workdays) .& :isolated)
    )
end

function n_tests(pop::Population, type::String)
    DataFrames.nrow(
        DataFramesMeta.@where(
            get_test_logs(pop.individuals),
            :type .== type
        )
    )
end

function n_infectious_per_day(pop::Population)
    tbl_status = get_status_over_time(pop)
    res = zeros(DataFrames.nrow(tbl_status))
    if "infectious, asymptomatic" in names(tbl_status)
        res += tbl_status[!, :"infectious, asymptomatic"]
    end
    if "infectious, symptomatic" in names(tbl_status)
        res += tbl_status[!, :"infectious, symptomatic"]
    end
    res
end

mean_contacts_per_day(pop::Population) = get_mean_contacts_over_time(pop)[!, :contacts]

function evaluate(pop::Population; tests::Vector{String} = ["pcr"], workdays::Vector{Int} = collect(0:4))
    n_infectious = n_infectious_per_day(pop)
    mean_contacts = mean_contacts_per_day(pop)
    # reduce to working days
    mean_contacts = mean_contacts[
        is_workday.(0:(length(mean_contacts) - 1); workdays = workdays)
    ]
    tbl = DataFrames.DataFrame(
        id = string(uuid4()),
        n_infected = n_infected(pop),
        workdays_missed = n_workdays_missed(pop; workdays = workdays),
        mean_daily_infectious = Statistics.mean(n_infectious),
        std_daily_infectious = Statistics.std(n_infectious),
        mean_daily_pp_contacts = Statistics.mean(mean_contacts)
    )
    for test in tests
        DataFrames.insertcols!(tbl, Symbol("n_$(test)_tests") => n_tests(pop, test))
    end
    return tbl
end

function evaluate(pops::Vector{T}; tests::Vector{String} = ["pcr"], workdays::Vector{Int} = collect(0:4)) where {T<:Population}
    vcat(evaluate.(pops; tests = tests, workdays = workdays)...)
end


function get_adjacency_matrix(pop::Population)
    n = length(pop.individuals)
    A = Matrix{Float64}(undef, n, n)
    for i = 1:n
        A[i, i] = Inf64
    end
    for i = 1:(n - 1), j = (i + 1):n
        a = pop.individuals[i]
        b = pop.individuals[j]
        A[i, j] = 0.0
        for g in pop.groups
            if is_member(a, g) & is_member(b, g)
                A[i, j] += g.meeting_probability
            end
        end
        A[j, i] = A[i, j]
    end
    return A
end
