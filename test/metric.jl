using Tullio
using Distributions
using LinearAlgebra
using Zygote
@testset "Trace of the metric" begin
    
    for n in 1:5

        #Get some coordiantes at random 
        r = rand(Uniform(3.0,1e5))      # Radial coordinate. 3.0 as rough lower limit of an event horizon
        θ = rand(Uniform(0.0, 2.0*π))   # Polar coord
        a = rand(Uniform(-0.99, 0.99))  # Spin parameter
        
        #Calculate the metric 
        coords = [0.0,r,θ,0.0]
        g = RelativisticDynamics.covariant_metric(coords,a)
        g_inverse = RelativisticDynamics.contravariant_metric(coords,a)
        
        #Check metric trace is OK 
        @tullio δ[a,c] := g[a,b] * g_inverse[b,c]  #:= allocates a new array
        @test isapprox(tr(δ),4.0)
    end

end



@testset "Christoffel tensor components" begin
    

    for n in 1:5

        #Get some coordiantes at random 
        t = rand(Uniform(3.0,1e5))       # Time coordinate - arbitratry since metric is time-independent
        r = rand(Uniform(3.0,1e5))       # Radial coordinate. 3.0 as rough lower limit of an event horizon
        θ = rand(Uniform(0.0, 2.0*π))    # Polar coordinate
        ϕ = rand(Uniform(0.0, 2.0*π))    # Azimuth coordinate - arbitratry since metric is axisymmetric
        a = rand(Uniform(-0.99, 0.99))   # Spin parameter

        coords = [t,r,θ,ϕ]
        
        #Calculate the metric 
        g = RelativisticDynamics.covariant_metric(coords,a)
        g_inverse = RelativisticDynamics.contravariant_metric(coords,a)

        #and the metric Jacobian using Zygote
        g_jacobian = jacobian(x -> RelativisticDynamics.covariant_metric_zygote(x,a), [t r θ ϕ])[1]
        g_jacobian = reshape(g_jacobian,(4,4,4)) # reshape

        #Christoffel
        @tullio Γ[μ,ν,λ] := 0.50*g_inverse[μ,ρ]*(g_jacobian[ρ,ν,λ] +g_jacobian[ρ,λ,ν]-g_jacobian[ν,λ,ρ])  #:= allocates a new array

        #Compare with the analytical solution
        Γ_analytical = RelativisticDynamics.christoffel(coords,a)
        @test isapprox(Γ,Γ_analytical)
    end
end 


@testset "Riemann tensor reduces to schwarzchild solution" begin



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

    for n in 1:5

        #Get some coordiantes at random 
        t = rand(Uniform(3.0,1e5))       # Time coordinate - arbitratry since metric is time-independent
        r = rand(Uniform(3.0,1e3))       # Radial coordinate. 3.0 as rough lower limit of an event horizon
        θ = rand(Uniform(0.0, 2.0*π))    # Polar coordinate
        ϕ = rand(Uniform(0.0, 2.0*π))    # Azimuth coordinate - arbitratry since metric is axisymmetric
        a = 0.0                          # Schwarzchild

 
        coords = [t,r,θ,ϕ]

        #Calculate the metric 
        g = RelativisticDynamics.covariant_metric(coords,a)
        g_inverse = RelativisticDynamics.contravariant_metric(coords,a)



        #Get the Riemann tensor solution, analytically
        Riemann_analytical = RelativisticDynamics.riemann(coords,a)
        

        #Get the schwazchild solution
        Riemann_covar_schwarzchild = schwarzchild_covariant_riemann(coords,a)

        #The contra/covar version of the Riemann tensor for schwarzchild
        @tullio Riemann_schwarzchild[μ,ν,ρ,σ] := g_inverse[μ,λ]*Riemann_covar_schwarzchild[λ,ν,ρ,σ]
       # end

        #Compare 
        @test isapprox(Riemann_analytical,Riemann_schwarzchild)

       
    end

end 

@testset "Kretschman scalar" begin
    

    """
    K =  Kretschmann_scalar(r,θ,a)
    Kretschman scalar for the Kerr metric
    """
    function Kretschmann_scalar(r,θ,a)
        Σ = r^2 + a^2 * cos(θ)^2
        return 48.0*(2.0*r^2-Σ)*(Σ^2-16.0*r^2*a^2*cos(θ)^2)/Σ^6
    end 


    for n in 1:5
    
        #Get some coordiantes at random 
        t = rand(Uniform(3.0,1e5))       # Time coordinate - arbitratry since metric is time-independent
        r = rand(Uniform(3.0,1e3))       # Radial coordinate. 3.0 as rough lower limit of an event horizon
        θ = rand(Uniform(0.0, 2.0*π))    # Polar coordinate
        ϕ = rand(Uniform(0.0, 2.0*π))    # Azimuth coordinate - arbitratry since metric is axisymmetric
        a = rand(Uniform(-0.99, 0.99))   # Spin parameter

        coords = [t,r,θ,ϕ]

        #Calculate the metric 
        g = RelativisticDynamics.covariant_metric(coords,a)
        g_inverse = RelativisticDynamics.contravariant_metric(coords,a)

       
        #Rieman tensor, mixed indices
        Riemann = RelativisticDynamics.riemann(coords,a)
      

        #Fully covariant form 
        @tullio Riemann_covar[μ,ν,ρ,σ] := g[μ,λ]*Riemann[λ,ν,ρ,σ]
        
        #Fully contravariant form 
        @tullio Riemann_contra[μ,ν,α,β] := g_inverse[μ,i]*g_inverse[ν,j]*g_inverse[α,k]*g_inverse[β,l]*Riemann_covar[i,j,k,l]
        

        #Kretschman scalar 
        @tullio Kretschman = Riemann_contra[α,β,μ,ν] * Riemann_covar[α,β,μ,ν] 
 

        Kretschman_formula = Kretschmann_scalar(r,θ,a)
        @test isapprox(Kretschman,Kretschman_formula)

    end

end 




















#EVERYTHING ABOVE THIS LINE WORKS--------------------------------


# @testset "Riemann tensor components" begin
    

#     #for n in 1:5
#     for n in 1:1

#         #Get some coordiantes at random 
#         t = rand(Uniform(3.0,1e5))       # Time coordinate - arbitratry since metric is time-independent
#         r = rand(Uniform(3.0,1e3))       # Radial coordinate. 3.0 as rough lower limit of an event horizon
#         θ = rand(Uniform(0.0, 2.0*π))    # Polar coordinate
#         ϕ = rand(Uniform(0.0, 2.0*π))    # Azimuth coordinate - arbitratry since metric is axisymmetric
#         a = rand(Uniform(-0.99, 0.99))   # Spin parameter

#         coords = [t,r,θ,ϕ]

#         #Calculate the metric 
#         println("metric")
#         g = RelativisticDynamics.covariant_metric(coords,a)
#         g_inverse = RelativisticDynamics.contravariant_metric(coords,a)

#         # First derivatives of the covariant metric 
#         println("metric derivative")
#         g_jacobian = jacobian(x -> RelativisticDynamics.covariant_metric_zygote(x,a), [t r θ ϕ])[1]
#         g_jacobian = reshape(g_jacobian,(4,4,4)) # reshape


#         # First derivatives of the contravariant metric 
#         println("inverse metric derivatieve")
#         g_inverse_jacobian = jacobian(x -> RelativisticDynamics.contravariant_metric(x,a), [t r θ ϕ])[1]
#         g_inverse_jacobian = reshape(g_inverse_jacobian,(4,4,4)) # reshape

#         display(g_inverse_jacobian)
#         # Second derivative of the metric 
#         println("hessian")
#         g_hessian = Zygote.hessian(x -> RelativisticDynamics.covariant_metric_zygote(x,a), [t r θ ϕ])[1]
#         g_hessian = reshape(g_hessian,(4,4,4,4)) # reshape

#         display(g_hessian)

#         # e.g. https://math.stackexchange.com/questions/628863/riemann-tensor-in-terms-of-the-metric-tensor
#         #This is R^{\rho}_{\sigma \mu \nu} 

#         @tullio Riemann[ρ,σ,μ,ν] := 
#         begin             
#         0.50*(g_inverse_∂[ρ,d,μ]*(g_∂[d,ν,σ] + g_∂[d,σ,ν] - g_∂[ν,σ,d])  + 
#         g_inverse[ρ,d]*(g_∂∂[d,ν,σ,μ]+g_∂∂[d,σ,ν,μ]-g_∂∂[ν,σ,d,μ]) -
#         g_inverse_∂[ρ,e,ν]*(g_∂[e,μ,σ] + g_∂[e,σ,μ] - g_∂[μ,σ,e])  -
#         g_inverse[ρ,e]*(g_∂∂[e,μ,σ,ν]+g_∂∂[e,σ,μ,ν]-g_∂∂[μ,σ,e,ν])
#        ) 
#         end 
        
#     #    @tullio begin
#     #         Riemann[ρ,σ,μ,ν] := 
#     #                             0.50*(g_inverse_∂[ρ,d,μ]*(g_∂[d,ν,σ] + g_∂[d,σ,ν] - g_∂[ν,σ,d])  + 
#     #                                   g_inverse[ρ,d]*(g_∂∂[d,ν,σ,μ]+g_∂∂[d,σ,ν,μ]-g_∂∂[ν,σ,d,μ]) -
#     #                                   g_inverse_∂[ρ,e,ν]*(g_∂[e,μ,σ] + g_∂[e,σ,μ] - g_∂[μ,σ,e])  -
#     #                                   g_inverse[ρ,e]*(g_∂∂[e,μ,σ,ν]+g_∂∂[e,σ,μ,ν]-g_∂∂[μ,σ,e,ν])
#     #                                  ) +
#     #                             0.25*g_inverse[ρ,f]*(g_∂[f,μ,λ] + g_∂[f,λ,μ] - g_∂[μ,λ,f])*g_inverse[λ,h]*(g_∂[h,ν,σ] + g_∂[h,σ,ν] - g_∂[ν,σ,h]) -
#     #                             0.25*g_inverse[ρ,i]*(g_∂[i,ν,λ] + g_∂[i,λ,ν] - g_∂[ν,λ,i])*g_inverse[λ,k]*(g_∂[k,μ,σ] + g_∂[k,σ,μ] - g_∂[μ,σ,k])
                                
#     #     end



#         #Compare with the analytical solution
#         Riemann_analytical = RelativisticDynamics.riemann(coords,a)
#         @test isapprox(Riemann,Riemann_analytical)


#         #Also double check the covariant form
#         @tullio Riemann_covar[μ,ν,ρ,σ] := g[μ,λ]*Riemann[λ,ν,ρ,σ]
        
#         #The purely covariant version of the Riemann tensor 
#         @tullio Riemann_covar_analytical[μ,ν,ρ,σ] := g[μ,λ]*Riemann_analytical[λ,ν,ρ,σ]
        

#         @test isapprox(Riemann_covar,Riemann_covar_analytical)


#     end

# end 























# # @testset "Riemann tensor components" begin
    

# #     #for n in 1:5
# #     for n in 1:1

# #         #Get some coordiantes at random 
# #         t = rand(Uniform(3.0,1e5))       # Time coordinate - arbitratry since metric is time-independent
# #         r = rand(Uniform(3.0,1e3))       # Radial coordinate. 3.0 as rough lower limit of an event horizon
# #         θ = rand(Uniform(0.0, 2.0*π))    # Polar coordinate
# #         ϕ = rand(Uniform(0.0, 2.0*π))    # Azimuth coordinate - arbitratry since metric is axisymmetric
# #         a = rand(Uniform(-0.99, 0.99))   # Spin parameter

# #         x = [t r θ ϕ]
        
# #         #Calculate the metric 
# #         g = RelativisticDynamics.covariant_metric(x,a)
# #         g_inverse = RelativisticDynamics.contravariant_metric(x,a)

# #         #Use automatic diff to get the first derivatives of the covariant metric
# #         g_∂ = zeros(Float64,4,4,4) # a derivative tensor
# #         shadow = onehot(x)
# #         for i = 1:length(shadow)
# #             g_∂[:,:,i] = autodiff(Forward, RelativisticDynamics.covariant_metric, DuplicatedNoNeed, Duplicated(x, shadow[i]), Const(a))[1]
# #         end 

# #         #...and first derivatives of the contravariant metric 
# #         g_inverse_∂ = zeros(Float64,4,4,4) # a derivative tensor
# #         shadow = onehot(x)
# #         for i = 1:length(shadow)
# #             g_inverse_∂[:,:,i] = autodiff(Forward, RelativisticDynamics.contravariant_metric, DuplicatedNoNeed, Duplicated(x, shadow[i]), Const(a))[1]
# #         end 

 
# #         # Second derivative of the metric 
# #         g_∂∂ = zeros(Float64,4,4,4,4) # a derivative tensor
        
# #         shadow = onehot(x)
# #         for i = 1:length(shadow)
# #             g_∂∂[:,:,:,i] = autodiff(Forward, g_∂, DuplicatedNoNeed, Duplicated(x, shadow[i]), Const(a))[1]
# #         end 






# #         # g_∂∂[1,1,:,:] = hessian(x -> RelativisticDynamics.metric_g11(x,a), [t r θ ϕ])
# #         # g_∂∂[2,2,:,:] = hessian(x -> RelativisticDynamics.metric_g22(x,a), [t r θ ϕ])
# #         # g_∂∂[3,3,:,:] = hessian(x -> RelativisticDynamics.metric_g33(x,a), [t r θ ϕ])
# #         # g_∂∂[4,4,:,:] = hessian(x -> RelativisticDynamics.metric_g44(x,a), [t r θ ϕ])
# #         # g_∂∂[1,4,:,:] = hessian(x -> RelativisticDynamics.metric_g14(x,a), [t r θ ϕ])
# #         # g_∂∂[4,1,:,:] = g_∂∂[1,4,:,:]


# #         # e.g. https://math.stackexchange.com/questions/628863/riemann-tensor-in-terms-of-the-metric-tensor
# #         #This is R^{\rho}_{\sigma \mu \nu} 
# #         @tensor begin
# #             Riemann[ρ,σ,μ,ν] := 
# #                                 0.50*(g_inverse_∂[ρ,d,μ]*(g_∂[d,ν,σ] + g_∂[d,σ,ν] - g_∂[ν,σ,d])  + 
# #                                       g_inverse[ρ,d]*(g_∂∂[d,ν,σ,μ]+g_∂∂[d,σ,ν,μ]-g_∂∂[ν,σ,d,μ]) -
# #                                       g_inverse_∂[ρ,e,ν]*(g_∂[e,μ,σ] + g_∂[e,σ,μ] - g_∂[μ,σ,e])  -
# #                                       g_inverse[ρ,e]*(g_∂∂[e,μ,σ,ν]+g_∂∂[e,σ,μ,ν]-g_∂∂[μ,σ,e,ν])
# #                                      ) +
# #                                 0.25*g_inverse[ρ,f]*(g_∂[f,μ,λ] + g_∂[f,λ,μ] - g_∂[μ,λ,f])*g_inverse[λ,h]*(g_∂[h,ν,σ] + g_∂[h,σ,ν] - g_∂[ν,σ,h]) -
# #                                 0.25*g_inverse[ρ,i]*(g_∂[i,ν,λ] + g_∂[i,λ,ν] - g_∂[ν,λ,i])*g_inverse[λ,k]*(g_∂[k,μ,σ] + g_∂[k,σ,μ] - g_∂[μ,σ,k])
                                
# #         end



# #         #Compare with the analytical solution
# #         Riemann_analytical = RelativisticDynamics.riemann(coords,a)
# #         @test isapprox(Riemann,Riemann_analytical)


# #         #Also double check the covariant form
# #         @tensor begin
# #             Riemann_covar[μ,ν,ρ,σ] := g[μ,λ]*Riemann[λ,ν,ρ,σ]
# #         end
# #         #The purely covariant version of the Riemann tensor 
# #         @tensor begin
# #             Riemann_covar_analytical[μ,ν,ρ,σ] := g[μ,λ]*Riemann_analytical[λ,ν,ρ,σ]
# #         end
# #         @test isapprox(Riemann_covar,Riemann_covar_analytical)


# #     end

# # end 



 












