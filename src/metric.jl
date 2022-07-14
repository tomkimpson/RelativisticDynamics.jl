function covariant_metric(r,θ,a)

    Δ = delta(r,a)
    Σ = sigma(r,θ,a)

    metric_covar = zeros(Float64,4,4)

    metric_covar[1,1] = 1.0 - 2.0*r / Σ
    metric_covar[2,2] = Σ / Δ
    metric_covar[3,3] = Σ
    metric_covar[4,4] = sin(θ)^2 * ((r^2 +a^2)^2 - Δ*a^2*sin(θ)^2) / Σ
    metric_covar[1,4] = -2.0*a*r*sin(θ)^2/Σ
    metric_covar[4,1] = metric_covar[1,4]
    return metric_covar
end 



function contravariant_metric(metric_covar,demoninator)

    

    metric_contra = zeros(Float64,4,4)
    

    metric_contra[1,1] = -metric_covar[4,4] / demoninator
    metric_contra[2,2] = 1.0 / metric_covar[2,2]
    metric_contra[3,3] = 1.0 / metric_covar[3,3]
    metric_contra[4,4] = -metric_covar[1,1] / demoninator

    metric_contra[1,4] = metric_covar[4,1]/demoninator
    metric_contra[4,1] = metric_contra[1,4]


    return metric_contra
end 