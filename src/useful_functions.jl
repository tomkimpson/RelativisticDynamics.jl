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
Convert a vector from contravariant form to convariant form using the metric 
"""
function convert_to_covariant(metric,vector)


    #vector_covar = vector
    #vector_covar[1] = metric[1,1]*vector[1] + metric[1,2]*vector[2] + metric[1,3]*vector[3] + metric[1,4]*vector[4]
    #vector_covar[2] = metric[2,1]*vector[1] + metric[2,2]*vector[2] + metric[2,3]*vector[3] + metric[2,4]*vector[4]
    #vector_covar[3] = metric[3,1]*vector[1] + metric[3,2]*vector[2] + metric[3,3]*vector[3] + metric[3,4]*vector[4]
    #vector_covar[4] = metric[4,1]*vector[1] + metric[4,2]*vector[2] + metric[4,3]*vector[3] + metric[4,4]*vector[4]

    @tensor begin
        vector_covar[μ] := metric[μ,ν] * vector[ν]  #:= allocates a new array
    end


    return vector_covar  

end 

"""
Kretschman scalar for the Kerr metric
"""
function Kretschmann_scalar(r,θ,a)

    Σ = sigma(r,θ,a)

    return 48.0*(2.0*r^2-Σ)*(Σ^2-16.0*r^2*a^2*cos(θ)^2)/Σ^6

end 


"""
Calculate the contravariant spin tensor
"""
function spintensor(xvector,pvector,svector,a,m0,ϵ)


    t,r,θ,ϕ = xvector
    Σ = sigma(r,θ,a)
    
    metric_trace =-sin(θ)^2*Σ^2

    #https://mathworld.wolfram.com/PermutationTensor.html
    permutation_tensor = ϵ/sqrt(abs(metric_trace))

    @tensor begin
         Stensor[μ,ν] := permutation_tensor[μ,ν,α,β]*pvector[α]*svector[β]/m0
    end

    return Stensor 

end 






