function covariant_metric(r,θ,a)

    Δ = delta(r,a)
    Σ = sigma(r,θ,a)

    metric_covar = zeros(Float64,4,4)

    metric_covar[1,1] = 1.0 - 2.0*r / Σ


    return metric_covar
end 



function contravariant_metric(r,θ,a)

    Δ = delta(r,a)
    Σ = sigma(r,θ,a)

    metric_contra = zeros(Float64,4,4)

    metric_contra[1,1] = 1.0 - 2.0*r / Σ


    return metric_contra
end 