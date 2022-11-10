using Zygote
using PreallocationTools


"""
    g=covariant_metric(coords,a)
Construct the NxN matrix of the covariant metric.
"""
function covariant_metric(coords,a)

    g = zeros(typeof(a),4,4)

    g[1,1] =   metric_g11(coords,a) 
    g[2,2] =   metric_g22(coords,a) 
    g[3,3] =   metric_g33(coords,a) 
    g[4,4] =   metric_g44(coords,a) 
    g[1,4] =   metric_g14(coords,a) 
    g[4,1] =   g[1,4] 

    return g

end 


"""
    g=contravariant_metric(coords,a)
Construct the NxN matrix of the contravariant metric.
Metric components are defined via indvidual functions to allow for auto diff in unit tests
"""
function contravariant_metric(coords,a)

    metric_contra = zeros(typeof(a),4,4)

    metric_contra[1,1] = metric_contra_g11(coords,a)
    metric_contra[2,2] = metric_contra_g22(coords,a)
    metric_contra[3,3] = metric_contra_g33(coords,a)
    metric_contra[4,4] = metric_contra_g44(coords,a)
    metric_contra[1,4] = metric_contra_g14(coords,a)
    metric_contra[4,1] = metric_contra[1,4]
    return metric_contra
end 




"""
    Γ = christoffel(coords,a)
The christoffel symbols of the Kerr metric. See e.g. https://arxiv.org/pdf/0904.4184.pdf
"""
function christoffel(coords,a)

    #See Catalogue of spacetimes: https://arxiv.org/pdf/0904.4184.pdf
    Γ = zeros(typeof(a),4,4,4) #Γ^{a}_{bc}
    r = coords[2]
    θ = coords[3]


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

    #set_to_zero!(Γ)
    return Γ
end 


"""
    R = riemann(coords,a) 
Riemann tensor components of the Kerr metric. First index is the contravariant, others are covariant   
"""
function riemann(coords,a)

    r = coords[2]
    θ = coords[3]
    Rtensor = zeros(typeof(a),4,4,4,4)

    Δ = delta(r,a)
    Σ = sigma(r,θ,a)

    #-----------------------------------------------------------------------#
    Rtensor[1,1,1,4] = 2.0*sin(θ)^2*a*r^2*(r^2-3.0*a^2*cos(θ)^2)/Σ^4
    Rtensor[1,1,2,3] = -2.0*(3.0*r^2-a^2*cos(θ)^2)*r*cos(θ)*a^2*sin(θ)/(Σ^3*Δ)
    Rtensor[1,2,1,2] = r*(2.0*(r^2+a^2)+a^2*sin(θ)^2)*(r^2-3.0*a^2*cos(θ)^2)/(Σ^3*Δ)
    Rtensor[1,2,1,3] = -(3.0*r^2-a^2*cos(θ)^2)*(3.0*(r^2+a^2)-2*r)*a^2*sin(θ)*cos(θ)/(Σ^3*Δ)
    Rtensor[1,2,2,4] = 3.0*sin(θ)^2*a*r*(r^2+a^2)*(r^2-3.0*a^2*cos(θ)^2)/(Σ^3*Δ)
    Rtensor[1,2,3,4] = -cos(θ)*sin(θ)*a*(3.0*r^2-a^2*cos(θ)^2)*(2.0*(r^2+a^2)^2+a^2*sin(θ)^2*Δ)/(Σ^3*Δ)
    
    #Antisymmetric under exchange of last two indices. 
    Rtensor[1,1,4,1] = -Rtensor[1,1,1,4]
    Rtensor[1,1,3,2] = -Rtensor[1,1,2,3]
    Rtensor[1,2,2,1] = -Rtensor[1,2,1,2]
    Rtensor[1,2,3,1] = -Rtensor[1,2,1,3] 
    Rtensor[1,2,4,2] = -Rtensor[1,2,2,4]  
    Rtensor[1,2,4,3] = -Rtensor[1,2,3,4]  


    #-----------------------------------------------------------------------#    
    

    Rtensor[1,3,1,2] = -(3.0*r^2-a^2*cos(θ)^2)*(3.0*(r^2+a^2)-4*r)*a^2*sin(θ)*cos(θ)/(Σ^3*Δ)
    Rtensor[1,3,1,3] = -r*(r^2-3.0*a^2*cos(θ)^2)*(r^2+a^2+2*a^2*sin(θ)^2)/Σ^3
    Rtensor[1,3,2,4] = -(3.0*r^2-a^2*cos(θ)^2)*((a^2+r^2)^2+2.0*a^2*sin(θ)^2*Δ)*a*sin(θ)*cos(θ)/(Σ^3*Δ)
    Rtensor[1,3,3,4] = -3.0*sin(θ)^2*r*a*(a^2+r^2)*(r^2-3.0*a^2*cos(θ)^2)/Σ^3
    Rtensor[1,4,1,4] = -sin(θ)^2*r*(r^2+3*a^2*sin(θ)^2-3*a^2)*((a^2+r^2)^2-a^2*sin(θ)^2*Δ)/Σ^4
    Rtensor[1,4,2,3] = (3.0*r^2-a^2*cos(θ)^2)*((a^2+r^2)^2-a^2*sin(θ)^2*Δ)*a*sin(θ)*cos(θ)/(Σ^3*Δ)
    
    #Symmetries
    Rtensor[1,3,2,1] = -Rtensor[1,3,1,2]
    Rtensor[1,3,3,1] = -Rtensor[1,3,1,3]
    Rtensor[1,3,4,2] = -Rtensor[1,3,2,4]
    Rtensor[1,3,4,3] = -Rtensor[1,3,3,4] 
    Rtensor[1,4,4,1] = -Rtensor[1,4,1,4]  
    Rtensor[1,4,3,2] = -Rtensor[1,4,2,3]  

    
    #-----------------------------------------------------------------------# 
    
    
    Rtensor[2,1,1,2] = r*(r^2-3.0*a^2*cos(θ)^2)*(a^2*sin(θ)^2+2.0*Δ)/Σ^4
    Rtensor[2,1,1,3] = -3.0*Δ*a^2*(3.0*r^2-a^2*cos(θ)^2)*sin(θ)*cos(θ)/Σ^4
    Rtensor[2,1,2,4] = r*a*(3.0*(r^2+a^2)-4.0*r)*(r^2-3.0*a^2*cos(θ)^2)*sin(θ)^2/Σ^4
    Rtensor[2,1,3,4] = -a*(2.0*(r^2+a^2)+a^2*sin(θ)^2)*(3.0*r^2-a^2*cos(θ)^2)*Δ*cos(θ)*sin(θ)/Σ^4
    Rtensor[2,3,1,4] = -cos(θ)*sin(θ)*a*(3.0*r^2-a^2*cos(θ)^2)*Δ/Σ^3
    Rtensor[2,3,2,3] = -r*(r^2-3.0*a^2*cos(θ)^2)/Σ^2
    Rtensor[2,4,1,2] = -sin(θ)^2*r*a*(r^2-3.0*a^2*cos(θ)^2)*(3.0*(r^2+a^2)-4.0*r)/Σ^4
    Rtensor[2,4,1,3] = a*(3.0*r^2-a^2*cos(θ)^2)*(r^2+a^2+2.0*a^2*sin(θ)^2)*Δ*cos(θ)*sin(θ)/Σ^4
    Rtensor[2,4,2,4] = -sin(θ)^2*r*(r^2-3.0*a^2*cos(θ)^2)*((r^2+a^2)^2+2.0*a^2*Δ*sin(θ)^2)/Σ^4
    Rtensor[2,4,3,4] = 3.0*a^2*(r^2+a^2)*(3.0*r^2-a^2*cos(θ)^2)*Δ*cos(θ)*sin(θ)^3/Σ^4 
    
    #Symmetries
    Rtensor[2,1,2,1] = -Rtensor[2,1,1,2]  
    Rtensor[2,1,3,1] = -Rtensor[2,1,1,3]  
    Rtensor[2,1,4,2] = -Rtensor[2,1,2,4] 
    Rtensor[2,1,4,3] = -Rtensor[2,1,3,4] 
    Rtensor[2,3,4,1] = -Rtensor[2,3,1,4]  
    Rtensor[2,3,3,2] = -Rtensor[2,3,2,3]  
    Rtensor[2,4,2,1] = -Rtensor[2,4,1,2]  
    Rtensor[2,4,3,1] = -Rtensor[2,4,1,3]  
    Rtensor[2,4,4,2] = -Rtensor[2,4,2,4]  
    Rtensor[2,4,4,3] = -Rtensor[2,4,3,4]   
    
    #-----------------------------------------------------------------------#
    
    
    Rtensor[3,1,1,2] = -3.0*a^2*(3.0*r^2-a^2*cos(θ)^2)*sin(θ)*cos(θ)/Σ^4
    Rtensor[3,1,1,3] = -r*(r^2-3.0*a^2*cos(θ)^2)*(2.0*a^2*sin(θ)^2+Δ)/Σ^4
    Rtensor[3,1,2,4] = -a*(3.0*r^2-a^2*cos(θ)^2)*(r^2+a^2+2.0*a^2*sin(θ)^2)*cos(θ)*sin(θ)/Σ^4
    Rtensor[3,1,3,4] = -r*a*(r^2-3.0*a^2*cos(θ)^2)*(3.0*(r^2+a^2)-2.0*r)*sin(θ)^2/Σ^4
    Rtensor[3,2,1,4] = a*(3.0*r^2-a^2*cos(θ)^2)*cos(θ)*sin(θ)/Σ^3
    Rtensor[3,2,2,3] = r*(r^2-3.0*a^2*cos(θ)^2)/(Δ*Σ^2)
    Rtensor[3,4,1,2] = a*(2.0*(r^2+a^2)+a^2*sin(θ)^2)*(3.0*r^2-a^2*cos(θ)^2)*cos(θ)*sin(θ)/Σ^4
    Rtensor[3,4,1,3] = r*a*(r^2-3.0*a^2*cos(θ)^2)*(3.0*(r^2+a^2)-2.0*r)*sin(θ)^2/Σ^4
    Rtensor[3,4,2,4] = 3.0*a^2*(r^2+a^2)*(3.0*r^2-a^2*cos(θ)^2)*cos(θ)*sin(θ)^3/Σ^4
    Rtensor[3,4,3,4] = r*(r^2-3.0*a^2*cos(θ)^2)*(2.0*(r^2+a^2)^2+a^2*Δ*sin(θ)^2)*sin(θ)^2/Σ^4 
    
    #Symmetries
    Rtensor[3,1,2,1] = -Rtensor[3,1,1,2]  
    Rtensor[3,1,3,1] = -Rtensor[3,1,1,3]  
    Rtensor[3,1,4,2] = -Rtensor[3,1,2,4] 
    Rtensor[3,1,4,3] = -Rtensor[3,1,3,4] 
    Rtensor[3,2,4,1] = -Rtensor[3,2,1,4]  
    Rtensor[3,2,3,2] = -Rtensor[3,2,2,3]  
    Rtensor[3,4,2,1] = -Rtensor[3,4,1,2]  
    Rtensor[3,4,3,1] = -Rtensor[3,4,1,3]  
    Rtensor[3,4,4,2] = -Rtensor[3,4,2,4]  
    Rtensor[3,4,4,3] = -Rtensor[3,4,3,4]   
    
    #-----------------------------------------------------------------------#
        
    Rtensor[4,1,2,3] = -(3.0*r^2-a^2*cos(θ)^2)*(a^2*sin(θ)^2-Δ)*a*(cos(θ)/sin(θ))/(Δ*Σ^3)
    Rtensor[4,2,1,2] = 3.0*r*a*(r^2-3.0*a^2*cos(θ)^2)/(Δ*Σ^3)
    Rtensor[4,2,1,3] = -a*(3.0*r^2-a^2*cos(θ)^2)*(2.0*a^2*sin(θ)^2+Δ)*(cos(θ)/sin(θ))/(Δ*Σ^3)
    Rtensor[4,2,2,4] = r*(r^2-3.0*a^2*cos(θ)^2)*(r^2+a^2+2.0*a^2*sin(θ)^2)/(Δ*Σ^3)
    Rtensor[4,2,3,4] = -(3.0*r^2-a^2*cos(θ)^2)*(3.0*(r^2+a^2)-2.0*r)*a^2*sin(θ)*cos(θ)/(Δ*Σ^3)
    Rtensor[4,3,1,2] = -(3.0*r^2-a^2*cos(θ)^2)*(a^2*sin(θ)^2+2.0*Δ)*a*(cos(θ)/sin(θ))/(Δ*Σ^3)
    Rtensor[4,3,1,3] = -3.0*r*a*(r^2-3.0*a^2*cos(θ)^2)/Σ^3
    Rtensor[4,3,2,4] = -a^2*(3*r^2-a^2*cos(θ)^2)*(3.0*(r^2+a^2)-4.0*r)*sin(θ)*cos(θ)/(Δ*Σ^3)
    Rtensor[4,3,3,4] = -r*(2.0*(r^2+a^2)+a^2*sin(θ)^2)*(r^2-3.0*a^2*cos(θ)^2)/Σ^3
    Rtensor[4,4,1,4] = -2.0*a*r^2*(r^2-3.0*a^2*cos(θ)^2)*sin(θ)^2/Σ^4
    Rtensor[4,4,2,3] = 2.0*r*a^2*(3.0*r^2-a^2*cos(θ)^2)*cos(θ)*sin(θ)/(Δ*Σ^3)
    
    #Symmetries
    Rtensor[4,1,3,2] = -Rtensor[4,1,2,3]  
    Rtensor[4,2,2,1] = -Rtensor[4,2,1,2]  
    Rtensor[4,2,3,1] = -Rtensor[4,2,1,3] 
    Rtensor[4,2,4,2] = -Rtensor[4,2,2,4] 
    Rtensor[4,2,4,3] = -Rtensor[4,2,3,4]  
    Rtensor[4,3,2,1] = -Rtensor[4,3,1,2]  
    Rtensor[4,3,3,1] = -Rtensor[4,3,1,3]  
    Rtensor[4,3,4,2] = -Rtensor[4,3,2,4]  
    Rtensor[4,3,4,3] = -Rtensor[4,3,3,4]  
    Rtensor[4,4,4,1] = -Rtensor[4,4,1,4]   
    Rtensor[4,4,3,2] = -Rtensor[4,4,2,3]   



    #Additional complicated terms - can these be algebrically simplified?
    part1 = -4.0*a^2*r^2*(a^2+r^2)*cos(θ)^2
    part2 = -4.0*a^2*(r^2-0.50*Σ)^2*sin(θ)^2 
    part3 = +2.0*a^2*r*cos(θ)^2*(Σ^2+2.0*a^2*r*sin(θ)^2)
    part4 = -(2.0*r^2-Σ)*(r*Σ^2-r^2*(a^2+r^2+Σ) +a^4*cos(θ)^2*sin(θ)^2)
    Rtensor[4,1,1,4]  =(part1+part2+part3+part4)/Σ^5
    Rtensor[4,1,4,1] = -Rtensor[4,1,1,4]
    

    #set_to_zero!(Rtensor)
    return Rtensor

end 


"""
    R = schwarzchild_covariant_riemann(coords,a)

Special case - the fully covariant components of the Riemann tensor for schwarzchild metric
Used for testing 
"""
function schwarzchild_covariant_riemann(coords,a)

    r = coords[2]
    θ = coords[3]
   
    Rtensor = zeros(typeof(a),4,4,4,4)

    Rtensor[1,2,1,2] = -2.0/r^3
    Rtensor[1,3,1,3] = (r-2.0)/r^2
    Rtensor[1,4,1,4] = (r-2.0)*sin(θ)^2/r^2
    Rtensor[2,3,2,3] = -1.0/(r-2.0)
    Rtensor[2,4,2,4] = -sin(θ)^2/(r-2.0)
    Rtensor[3,4,3,4] = r*2.0*sin(θ)^2


    #Symmetries
    Rtensor[1,2,2,1] = -Rtensor[1,2,1,2]
    Rtensor[1,3,3,1] = -Rtensor[1,3,1,3]
    Rtensor[1,4,4,1] = -Rtensor[1,4,1,4]
    Rtensor[2,3,3,2]=  -Rtensor[2,3,2,3]
    Rtensor[2,4,4,2] = -Rtensor[2,4,2,4]
    Rtensor[3,4,4,3] = -Rtensor[3,4,3,4]

    Rtensor[2,1,1,2] = -Rtensor[1,2,1,2]
    Rtensor[3,1,1,3] = -Rtensor[1,3,1,3]
    Rtensor[4,1,1,4] = -Rtensor[1,4,1,4]
    Rtensor[3,2,2,3]=  -Rtensor[2,3,2,3]
    Rtensor[4,2,2,4] = -Rtensor[2,4,2,4]
    Rtensor[4,3,3,4] = -Rtensor[3,4,3,4]

    Rtensor[2,1,2,1] = Rtensor[1,2,1,2]
    Rtensor[3,1,3,1] = Rtensor[1,3,1,3] 
    Rtensor[4,1,4,1] = Rtensor[1,4,1,4] 
    Rtensor[3,2,3,2] = Rtensor[2,3,2,3] 
    Rtensor[4,2,4,2] = Rtensor[2,4,2,4] 
    Rtensor[4,3,4,3] = Rtensor[3,4,3,4] 

    return Rtensor



end 

"""
    g = metric_g11(coords,a)
The covariant tt component of the Kerr metric in Boyer Lindquist coordinates
"""
function metric_g11(x,a)
    t,r,θ,ϕ =  x[1],x[2],x[3],x[4]
    Σ = RelativisticDynamics.sigma(r,θ,a)
    return -(1.0 - 2.0*r / Σ)
end 


"""
    g = metric_g22(coords,a)
The covariant rr component of the Kerr metric in Boyer Lindquist coordinates
"""
function metric_g22(x,a)
    t,r,θ,ϕ =  x[1],x[2],x[3],x[4]
    Σ = RelativisticDynamics.sigma(r,θ,a)
    Δ = RelativisticDynamics.delta(r,a)
    return Σ / Δ
end 

"""
    g = metric_g33(coords,a)
The covariant θθ component of the Kerr metric in Boyer Lindquist coordinates
"""
function metric_g33(x,a)
    t,r,θ,ϕ =  x[1],x[2],x[3],x[4]
    Σ = RelativisticDynamics.sigma(r,θ,a)
    return Σ 
end 

"""
    g = metric_g44(coords,a)
The covariant ϕϕ component of the Kerr metric in Boyer Lindquist coordinates
"""
function metric_g44(x,a)
    t,r,θ,ϕ =  x[1],x[2],x[3],x[4]
    Σ = RelativisticDynamics.sigma(r,θ,a)
    Δ = RelativisticDynamics.delta(r,a)
    return sin(θ)^2 * ((r^2 +a^2)^2 - Δ*a^2*sin(θ)^2) / Σ
end 

"""
    g = metric_g14(coords,a)
The covariant tϕ component of the Kerr metric in Boyer Lindquist coordinates
"""
function metric_g14(x,a)
    t,r,θ,ϕ =  x[1],x[2],x[3],x[4]
    Σ = RelativisticDynamics.sigma(r,θ,a)
    return -2.0*a*r*sin(θ)^2/Σ
end 


"""
    g = metric_g11(coords,a)
The covariant tt component of the Kerr metric in Boyer Lindquist coordinates
"""
function metric_contra_g11(x,a)
    t,r,θ,ϕ =  x[1],x[2],x[3],x[4]
    Σ = RelativisticDynamics.sigma(r,θ,a)
    Δ = RelativisticDynamics.delta(r,a)
    covar = sin(θ)^2 * ((r^2 +a^2)^2 - Δ*a^2*sin(θ)^2) / Σ
    denom = Δ*sin(θ)^2
    return -covar/denom 
end 


"""
    g = metric_g22(coords,a)
The covariant rr component of the Kerr metric in Boyer Lindquist coordinates
"""
function metric_contra_g22(x,a)
    t,r,θ,ϕ =  x[1],x[2],x[3],x[4]
    Σ = RelativisticDynamics.sigma(r,θ,a)
    Δ = RelativisticDynamics.delta(r,a)
    return  Δ/Σ 
end 

"""
    g = metric_g33(coords,a)
The covariant θθ component of the Kerr metric in Boyer Lindquist coordinates
"""
function metric_contra_g33(x,a)
    t,r,θ,ϕ =  x[1],x[2],x[3],x[4]
    Σ = RelativisticDynamics.sigma(r,θ,a)
    return 1.0/Σ 
end 

"""
    g = metric_g44(coords,a)
The covariant ϕϕ component of the Kerr metric in Boyer Lindquist coordinates
"""
function metric_contra_g44(x,a)
    t,r,θ,ϕ =  x[1],x[2],x[3],x[4]
    Σ = RelativisticDynamics.sigma(r,θ,a)
    Δ = RelativisticDynamics.delta(r,a)
    covar = -(1.0 - 2.0*r / Σ)
    denom = Δ*sin(θ)^2

    return -covar/denom
end 

"""
    g = metric_g14(coords,a)
The covariant tϕ component of the Kerr metric in Boyer Lindquist coordinates
"""
function metric_contra_g14(x,a)
    t,r,θ,ϕ =  x[1],x[2],x[3],x[4]
    Σ = RelativisticDynamics.sigma(r,θ,a)
    Δ = RelativisticDynamics.delta(r,a)
    covar = -2.0*a*r*sin(θ)^2/Σ
    denom = Δ*sin(θ)^2
    return covar/denom
end 







