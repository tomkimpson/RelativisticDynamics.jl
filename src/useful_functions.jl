using LinearAlgebra
using ChainRulesCore
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
    p_{μ} = convert_to_covariant(metric,p^{μ})
Convert a vector from contravariant form to convariant form using the covariant metric 
"""
function convert_to_covariant(metric,vector)

    @tullio vector_covar[μ] := metric[μ,ν] * vector[ν]  #:= allocates a new array
 
    return vector_covar  

end 



"""
    K =  Kretschmann_scalar(r,θ,a)
Kretschman scalar for the Kerr metric
"""
function Kretschmann_scalar(r,θ,a)

    Σ = sigma(r,θ,a)

    return 48.0*(2.0*r^2-Σ)*(Σ^2-16.0*r^2*a^2*cos(θ)^2)/Σ^6

end 


"""
    l = calculate_levi(NF)
Determine the Levi-civita psuedo tensor
"""
function calculate_levi()

    # Levi civita tensor
    levi = zeros(Integer,4,4,4,4) 
    ChainRulesCore.ignore_derivatives() do # This can be safely ignored by the differentiator - no dependence on the input parameters.

        for i in 1:4
            for j in 1:4
                for k in 1:4
                    for l in 1:4
                        permutation_vector = [i,j,k,l]
                        levi[i,j,k,l] = levicivita(permutation_vector) #This is [i,j,k,l] from e.g. https://mathworld.wolfram.com/PermutationTensor.html
                    end
                end 
            end
        end 

    end  

    return levi 
end 





"""
    ϵ = permutation_tensor(metric)
Calcualte the Levi-civita tensor in an arbitrary basis.
"""
function permutation_tensor(metric)

    ϵ = calculate_levi()

    det_g = det(metric)

    return ϵ/sqrt(abs(det_g))

end 


"""
    S^{μ ν} = spintensor(levi,pvector,svector,m0)
Calculate the contravariant spin tensor. 
"""
function spintensor(levi,pvector,svector,m0)

    @tullio Stensor[μ,ν] := levi[μ,ν,α,β]*pvector[α]*svector[β]/m0
    
    return Stensor 

end 






