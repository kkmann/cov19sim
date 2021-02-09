struct DoNothing <: Policy end

test_and_isolate!(pol::DoNothing, gr::Group) = return
