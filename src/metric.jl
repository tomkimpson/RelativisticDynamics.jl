"""
    Δ = delta(r,a)
The well-known delta function of the Kerr metric
"""
function delta(r,a)
return r^2 -2.0*r + a^2
end 

"""
    Σ = sigma(r,θ,a)
The well-known sigma function of the Kerr metric
"""
function sigma(r,θ,a)
return r^2 + a^2 * cos(θ)^2
end 


# """
#     g=covariant_metric(coords,a)
# Construct the NxN matrix of the covariant metric.
# """
# function covariant_metric(coords,a)

#     #xs = zeros(typeof(a),4,4)
#     #g = Zygote.Buffer(xs) #zeros(typeof(a),4,4)
    
#     g = zeros(typeof(a),4,4)

#     t,r,θ,ϕ =  coords[1],coords[2],coords[3],coords[4]
#     Σ = sigma(r,θ,a)
#     Δ = delta(r,a)


#     g[1,1] =   -(1.0 - 2.0*r / Σ)
#     g[2,2] =   Σ / Δ
#     g[3,3] =   Σ 
#     g[4,4] =   sin(θ)^2 * ((r^2 +a^2)^2 - Δ*a^2*sin(θ)^2) / Σ
#     g[1,4] =   -2.0*a*r*sin(θ)^2/Σ 
#     g[4,1] =   g[1,4] 

#    # println("Normal metric")
#    # display(g)

#     #return copy(g)
#     return g

# end 



function covariant_metric(coords,a)

    xs = zeros(typeof(a),4,4)
    g = Zygote.bufferfrom(xs) 

    t,r,θ,ϕ =  coords[1],coords[2],coords[3],coords[4]
    Σ = sigma(r,θ,a)
    Δ = delta(r,a)

    g[1,1] =   -(1.0 - 2.0*r / Σ)
    g[2,2] =   Σ / Δ
    g[3,3] =   Σ 
    g[4,4] =   sin(θ)^2 * ((r^2 +a^2)^2 - Δ*a^2*sin(θ)^2) / Σ
    g[1,4] =   -2.0*a*r*sin(θ)^2/Σ 
    g[4,1] =   g[1,4] 

    return copy(g)

end 



"""
    g=contravariant_metric(coords,a)
Construct the NxN matrix of the contravariant metric.
Metric components are defined via indvidual functions to allow for auto diff in unit tests
"""
function contravariant_metric(coords,a)



    xs = zeros(typeof(a),4,4)
    metric_contra = Zygote.bufferfrom(xs)
    t,r,θ,ϕ =  coords[1],coords[2],coords[3],coords[4]
    Σ = sigma(r,θ,a)
    Δ = delta(r,a)

    
    denom = Δ*sin(θ)^2

    metric_contra[1,1] = -(sin(θ)^2 * ((r^2 +a^2)^2 - Δ*a^2*sin(θ)^2) / Σ)/denom 
    metric_contra[2,2] = Δ/Σ
    metric_contra[3,3] = 1.0/Σ 
    metric_contra[4,4] = -(-(1.0 - 2.0*r / Σ))/denom
    metric_contra[1,4] = (-2.0*a*r*sin(θ)^2/Σ)/denom
    metric_contra[4,1] = metric_contra[1,4]
    return copy(metric_contra)
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
    

    
    return Rtensor

end 




"""
    E,L,Q = ELQ(a,α,e,ι,D)
Calculate the energy, angular momentum and Carter constant given the Keplerian orbital parameters, and the BH spin/direction
Specific to the Kerr metric, see Schmidt 2002 arXiv:0202090
"""
function ELQ(a,α,e,ι,D)



    #First define some nested functions
    function mapping_f(r,a,zminus)
        Δ = delta(r,a)
        return r^4 +a^2 * (r*(r+2.0) + zminus^2 * Δ)
    end

    function mapping_g(r,a)
        return 2*a*r
    end

    function mapping_h(r,a,zminus)
        Δ = delta(r,a)
        return r*(r-2.0) + (zminus^2 * Δ)/(1.0 - zminus^2)
    end

    function mapping_d(r,a,zminus)
        Δ = delta(r,a)
        return (r^2 +a^2 * zminus^2)*Δ
    end
        

    zm = cos(ι)

    r_periapsis= α*(1-e)
    r_apoapsis = α*(1+e)
    
    #Define some orbital coefficients
    f1 = mapping_f(r_periapsis,a,zm)
    g1 = mapping_g(r_periapsis,a)
    h1 = mapping_h(r_periapsis,a,zm)
    d1 = mapping_d(r_periapsis,a,zm)

    f2 = mapping_f(r_apoapsis,a,zm)
    g2 = mapping_g(r_apoapsis,a)
    h2 = mapping_h(r_apoapsis,a,zm)
    d2 = mapping_d(r_apoapsis,a,zm)

 
    #Determinants 
    κ = d1*h2 - d2*h1 
    ϵ = d1*g2 - d2*g1 
    ρ = f1*h2 - f2*h1
    η = f1*g2 - f2*g1 
    σ = g1*h2 - g2*h1 


    #Energy
    E_numerator   = κ*ρ+2.0*η*σ-2.0*D*sqrt(σ*(σ*ϵ^2 + ρ*ϵ*κ-η*κ^2))
    E_denominator = ρ^2 + 4.0*η*σ
    E = sqrt(E_numerator/E_denominator)

    #Angular momentum 
    L = -g1*E/h1 + D*sqrt(g1^2 * E^2 +(f1*E^2 -d1)*h1)/h1

    #Carter 
    Q = zm^2*(a^2*(1.0-E^2)+L^2/(1.0-zm^2))

    return E,L,Q


end 







