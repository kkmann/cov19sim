struct HeavyTailsModel{T,DM} <: DiseaseModel
    base_dm::DM # underlying basic disease model
    l::T # lengthscale of covariance function
    scale::T # standard deviation of covariance function
    df::T # degrees of freedom
end

# constructor with default values (except infection model)
function HeavyTailsModel(
    dm::DM;
    l::T = 3.0,
    scale::T = 1.0,
    df::T = 3.0
) where {T<:Real, DM<:DiseaseModel}
    HeavyTailsModel{T,DM}(dm, l, scale, df)
end

function sample(dm::HeavyTailsModel; n_retry = 1000)

    # sample mean from underlying model
    dt = sample(dm.base_dm)
    l10vl = log10.(dt.vl)
    if length(l10vl) < 3
        throw(InexactError)
    end

    # construct covariance matrix
    K = Matrix{eltype(l10vl)}(undef, length(l10vl), length(l10vl))
    for i = 1:length(l10vl)
        K[i, i] = dm.scale^2 + sqrt(eps())
        for j = 1:(i - 1)
            d = abs(i - j)
            K[i, j] = dm.scale^2*exp(-(d/dm.l)^2)
            K[j, i] = K[i, j]
        end
    end

    # condition on initial and last value being 0
    # see https://arxiv.org/abs/1402.4306
    idx_known = [1, length(l10vl)]
    idx_unknown = 2:(length(l10vl) - 1) # indices in reordered values
    K11 = K[idx_known, idx_known]
    K12 = K[idx_known, idx_unknown]
    K21 = K[idx_unknown, idx_known]
    K22 = K[idx_unknown, idx_unknown]
    KK_ = (dm.df - 2)/(dm.df + 2 - 2) * (K22 - K21*inv(Symmetric(K11))*K12)

    # sample from multivariate t distribution conditional on log10 VL + noise being positive
    X1 = Distributions.MvNormal(zeros(length(idx_unknown)), Symmetric(KK_))
    X2 = Distributions.Chisq(dm.df)
    tmp = zeros(length(l10vl))
    success = false
    for i = 1:n_retry
        tmp = l10vl + vcat(0.0, rand(X1, 1)[:, 1] ./ sqrt(rand(X2, 1)[1]/dm.df), 0.0)
        # check if we have sensible vl values
        if all(tmp .>= 0)
            success = true
            break
        end
    end
    if !success
        throw(InexactError)
    end

    dt.vl .= 10 .^ tmp
    return dt
end
