n = 5

gr = Group(
    [Individual(LarremoreModel(0.05), 0.01) for i = 1:5],
    DoNothing(),
    .33,
    collect(0:4)
)

A = get_adjacency_matrix(gr)
for i in 1:(n - 1), j = (i + 1):n
    @test A[i, j] == .33
end
@test all([A[i, i] for i = 1:n] .== Inf64)
