"""
    p_{μ} = convert_to_covariant(metric,p^{μ})
Convert a vector from contravariant form to convariant form using the covariant metric 
"""
function convert_to_covariant(metric,vector)

    @tullio vector_covar[μ] := metric[μ,ν] * vector[ν]  #:= allocates a new array
 
    return vector_covar  

end 


"""
    l = levi_civita_symbol(NF)
Determine the Levi-civita psuedo tensor
"""
function levi_civita_symbol()

    # Levi civita tensor
    levi = zeros(Integer,4,4,4,4) 
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

    return levi 
end 

"""
    ϵ = levi_civita_tensor(metric)
Calcualte the Levi-civita tensor in an arbitrary basis.
"""
function levi_civita_tensor(metric)

    ϵ = levi_civita_symbol()

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




"""
    CartesianCoordinates()
A struct to hold...
"""
@with_kw struct CartesianTrajectory           
    x  
    y  
    z  
end


"""
Given Boyer-Lindquist coordinates, return an struct of the Cartesian coordinates
"""
function boyer_lindquist_to_cartesian(solution,a)

    interpolation_factor = 10 #By default upsample/interpolate by a factor of 10 for smooth plotting 

    T = range(first(solution.t),last(solution.t),length=length(solution.t)*interpolation_factor)
    p = solution(T)

    # Extract relevant data from the interpolated solution 
    r = p[2,:]
    θ = p[3,:] 
    ϕ = p[4,:]

    # Boyer lindquist to Cartesian 
    w = sqrt.(r.^2 .+ a^2) 
    x = w .* sin.(θ) .* cos.(ϕ)
    y = w .* sin.(θ) .* sin.(ϕ)
    z = r .* cos.(θ)

   return CartesianTrajectory(x,y,z)

end 






