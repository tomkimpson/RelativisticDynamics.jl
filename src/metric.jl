function covariant_metric(r,θ,a)

    Δ = delta(r,a)
    Σ = sigma(r,θ,a)

    metric_covar = zeros(Float64,4,4)

    metric_covar[1,1] =  metric_gtt(r,θ,a) #-(1.0 - 2.0*r / Σ)
    metric_covar[2,2] =  metric_grr(r,θ,a) #Σ / Δ
    metric_covar[3,3] =  metric_gθθ(r,θ,a)
    metric_covar[4,4] =  metric_gϕϕ(r,θ,a) #sin(θ)^2 * ((r^2 +a^2)^2 - Δ*a^2*sin(θ)^2) / Σ
    metric_covar[1,4] =  metric_gtϕ(r,θ,a) #-2.0*a*r*sin(θ)^2/Σ
    metric_covar[4,1] =  metric_covar[1,4]
    return metric_covar
end 



function contravariant_metric(metric_covar,demoninator)

    #see https://www.roma1.infn.it/teongrav/onde19_20/kerr.pdf

    metric_contra = zeros(Float64,4,4)
    

    metric_contra[1,1] = -metric_covar[4,4] / demoninator
    metric_contra[2,2] = 1.0 / metric_covar[2,2]
    metric_contra[3,3] = 1.0 / metric_covar[3,3]
    metric_contra[4,4] = -metric_covar[1,1] / demoninator

    metric_contra[1,4] = metric_covar[4,1]/demoninator
    metric_contra[4,1] = metric_contra[1,4]


    return metric_contra
end 




"""
Christoffel symbols of the second kind 
"""
function christoffel(r,θ,a)

    #See Catalogue of spacetimes: https://arxiv.org/pdf/0904.4184.pdf
    Γ = zeros(Float64,4,4,4) #Γ^{a}_{bc}

    Δ = delta(r,a)
    Σ = sigma(r,θ,a)
    r2 = r^2 
    a2 = a^2
    c2 = cos(θ)^2
    s2 = sin(θ)^2
    
    
    #t
    Γ[1,1,2] = (r2 + a2)*(r2-a2*c2) / (Σ^2 * Δ)
    Γ[1,1,3] = -2.0*r*a2*sin(θ)*cos(θ)/Σ^2
    Γ[1,2,4] = -a*((3*r2-a2)*(r2+a2) -a2*(r2-a2)*s2)*s2/(Σ^2 * Δ)
    Γ[1,3,4] = 2.0*r*a^3*sin(θ)^3*cos(θ)/Σ^2

    #r
    Γ[2,1,1] = Δ*(r2 - a2*c2)/Σ^3
    Γ[2,1,4] = -a*Δ*(r2-a2*c2)*s2/Σ^3
    Γ[2,2,2] = (-(r2 -a2)+a2*s2*(r-1.0))/(Σ*Δ)
    Γ[2,2,3] = -a2*sin(θ)*cos(θ)/Σ
    Γ[2,3,3] = -r*Δ/Σ
    Γ[2,4,4] = Δ*s2*(-2.0*r*Σ^2 + 2.0*a2*s2*(r2-a2*c2))/(2.0*Σ^3)


    #θ
    Γ[3,1,1] = -2.0*r*a2*sin(θ)*cos(θ)/Σ^3
    Γ[3,1,4] = 2.0*r*a*(r2+a2)*sin(θ)*cos(θ)/Σ^3
    Γ[3,2,2] = a2*sin(θ)*cos(θ)/(Σ*Δ)
    Γ[3,2,3] = r/Σ
    Γ[3,3,3] = -a2*sin(θ)*cos(θ)/Σ
    Γ[3,4,4] = -sin(θ)*cos(θ)*(((r2+a2)^2 - a2*Δ*s2)*Σ + (r2+a2)*2.0*a2*r*s2)/Σ^3


    #Φ
    Γ[4,1,2] = a*(r2-a2*cos(θ)^2)/(Σ^2*Δ)
    Γ[4,1,3] = -2.0*r*a*(cos(θ)/sin(θ))/Σ^2
    Γ[4,2,4] = (2.0*r*Σ^2 + 2.0*(a2^2*s2*c2 - r2*(Σ+r2+a2)))/(2.0*Σ^2*Δ)
    Γ[4,3,4] = cot(θ)*(Σ^2 + 2.0*a2*r*s2)/Σ^2


    #By symmetry. https://mathworld.wolfram.com/ChristoffelSymboloftheSecondKind.html
    Γ[1,2,1] = Γ[1,1,2]
    Γ[1,3,1] = Γ[1,1,3]
    Γ[1,4,2] = Γ[1,2,4]
    Γ[1,4,3] = Γ[1,3,4]


    Γ[2,4,1] = Γ[2,1,4]
    Γ[2,3,2]  = Γ[2,2,3]

    Γ[3,4,1] = Γ[3,1,4] 
    Γ[3,3,2] = Γ[3,2,3] 

    Γ[4,2,1]=Γ[4,1,2]  
    Γ[4,3,1]=Γ[4,1,3]  
    Γ[4,4,2]=Γ[4,2,4]  
    Γ[4,4,3]=Γ[4,3,4]  


    return Γ
end 



# Pure function definitions
# Function definitions which do not mutate inputs.
# Allows for checks via AutoDif https://fluxml.ai/Zygote.jl/latest/limitations/
function metric_gtt(r,θ,a)
    Σ = sigma(r,θ,a)
    return -(1.0 - 2.0*r / Σ)
end 

function metric_grr(r,θ,a)
    Σ = sigma(r,θ,a)
    Δ = delta(r,a)
    return Σ / Δ
end 

function metric_gθθ(r,θ,a)
    Σ = sigma(r,θ,a)
    return Σ 
end 

function metric_gϕϕ(r,θ,a)
    Σ = sigma(r,θ,a)
    Δ = delta(r,a)
    return sin(θ)^2 * ((r^2 +a^2)^2 - Δ*a^2*sin(θ)^2) / Σ
end 

function metric_gtϕ(r,θ,a)
    Σ = sigma(r,θ,a)
    return -2.0*a*r*sin(θ)^2/Σ
end 
