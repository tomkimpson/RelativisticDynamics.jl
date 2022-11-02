using LinearAlgebra

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
Convert a vector from contravariant form to convariant form using the covariant metric 
"""
function convert_to_covariant(metric,vector)


    @tullio vector_covar[μ] := metric[μ,ν] * vector[ν]  #:= allocates a new array
 
    return vector_covar  

end 



"""
Kretschman scalar for the Kerr metric
"""
function Kretschmann_scalar(r,θ,a)

    Σ = sigma(r,θ,a)

    return 48.0*(2.0*r^2-Σ)*(Σ^2-16.0*r^2*a^2*cos(θ)^2)/Σ^6

end 



function permutation_tensor(metric,ϵ)

    det_g = det(metric)

    return ϵ/sqrt(abs(det_g))

end 


"""
Calculate the contravariant spin tensor. See e.g. https://mathworld.wolfram.com/PermutationTensor.html
"""
function spintensor(levi,pvector,svector,m0)


    @tullio Stensor[μ,ν] := levi[μ,ν,α,β]*pvector[α]*svector[β]/m0
    
    return Stensor 

end 






