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

"""
f = mapping_f(r,a,zminus)
Mapping function `f` used when converting from Keplerian orbital parameters to constants of motion
"""
function mapping_f(r,a,zminus)
Δ = delta(r,a)
return r^4 +a^2 * (r*(r+2.0) + zminus^2 * Δ)
end


"""
g = mapping_g(r,a)
Mapping function `g` used when converting from Keplerian orbital parameters to constants of motion
"""
function mapping_g(r,a)
return 2*a*r
end

"""
h = mapping_h(r,a,zminus)
Mapping function `h` used when converting from Keplerian orbital parameters to constants of motion
"""
function mapping_h(r,a,zminus)
Δ = delta(r,a)
return r*(r-2.0) + (zminus^2)/(1.0 - zminus^2) *Δ
end


"""
d = mapping_d(r,a,zminus)
Mapping function `d` used when converting from Keplerian orbital parameters to constants of motion
"""
function mapping_d(r,a,zminus)
Δ = delta(r,a)
return (r^2 +a^2 * zminus^2)*Δ
end



"""
Convert a vector from contravariant form to convariant form using the metric 
"""
function convert_to_covariant(metric,vector)


    vector_covar = vector
    vector_covar[1] = metric[1,1]*vector[1] + metric[1,2]*vector[2] + metric[1,3]*vector[3] + metric[1,4]*vector[4]
    vector_covar[2] = metric[2,1]*vector[1] + metric[2,2]*vector[2] + metric[2,3]*vector[3] + metric[2,4]*vector[4]
    vector_covar[3] = metric[3,1]*vector[1] + metric[3,2]*vector[2] + metric[3,3]*vector[3] + metric[3,4]*vector[4]
    vector_covar[4] = metric[4,1]*vector[1] + metric[4,2]*vector[2] + metric[4,3]*vector[3] + metric[4,4]*vector[4]

    return vector_covar  

end 


using Plots 
function BoyerLindquistPlot(solution,a)

    r = solution[2,:]
    θ = solution[3,:]
    ϕ = solution[4,:]


    w = sqrt.(r.^2 .+ a^2) 
    x = w .* sin.(θ) .* cos.(ϕ)
    y = w .* sin.(θ) .* sin.(ϕ)
    z = r .* cos.(θ)

    #Set the backend 
    #gr()
    
    #plotly(lw=3)
    #plot(x,y)
    scatter(x,y,z)
    #plot(solution,vars=[1])

end 