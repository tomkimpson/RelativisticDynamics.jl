@testset "Trace of the metric" begin
    

    for n in 1:3

        #Get some coordiantes at random 
        r = rand(Uniform(3.0,1e5))      # Radial coordinate. 3.0 as rough lower limit of an event horizon
        θ = rand(Uniform(0.0, 2.0*π))   # Polar coord
        a = rand(Uniform(-0.99, 0.99))  # Spin parameter
        
        #Calculate the metric 
        g = RelativisticDynamics.covariant_metric(r,θ,a)
        Δ = RelativisticDynamics.delta(r,a)
        g_inverse = RelativisticDynamics.contravariant_metric(g,Δ*sin(θ)^2)
        #Check metric trace is OK 
        @tensor begin
        δ[a,c] := g[a,b] * g_inverse[b,c]  #:= allocates a new array
        end
        @test tr(δ)==4.0 
    end


    # a=1
    # @test a==1


    # blob = gradient(x -> 3x^2 + 2x + 1, 5)
   
    # println(tr(blob))

    # Γ = RelativisticDynamics.christoffel(r,θ,a)
    
    #println(Γ)
end













@testset "Christoffel symbols" begin
    

    for n in 1:1

        #Get some coordiantes at random 
        #r = rand(Uniform(3.0,1e5))      # Radial coordinate. 3.0 as rough lower limit of an event horizon
        #θ = rand(Uniform(0.0, 2.0*π))   # Polar coord
        #a = rand(Uniform(-0.99, 0.99))  # Spin parameter


        t = 0.0
        r = 20.0
        θ = π/4.0
        ϕ = π/6.0 
        a = 0.0

        println(r," ",θ," ", a)

        #Calculate the metric 
        g = RelativisticDynamics.covariant_metric(r,θ,a)
        Δ = RelativisticDynamics.delta(r,a)
        g_inverse = RelativisticDynamics.contravariant_metric(g,Δ*sin(θ)^2)

        #Use automatic diff to get the first derivatives of the metric  of the metric 
        #Only calculate r and θ derivatives
        #`jacobian()` returns a tuple so [1] extracts the vector of floats 
        g_∂ = zeros(Float64,4,4,4) # a derivative tensor 

        g_∂[1,1,:] = jacobian(x -> RelativisticDynamics.metric_g11(x,a), [t r θ ϕ])[1]
        g_∂[2,2,:] = jacobian(x -> RelativisticDynamics.metric_g22(x,a), [t r θ ϕ])[1]
        g_∂[3,3,:] = jacobian(x -> RelativisticDynamics.metric_g33(x,a), [t r θ ϕ])[1]
        g_∂[4,4,:] = jacobian(x -> RelativisticDynamics.metric_g44(x,a), [t r θ ϕ])[1]
        g_∂[1,4,:] = jacobian(x -> RelativisticDynamics.metric_g14(x,a), [t r θ ϕ])[1]

        @tensor begin
            Γ[μ,ν,λ] := 0.50*g_inverse[μ,ρ]*(g_∂[ρ,ν,λ] +g_∂[ρ,λ,ν]-g_∂[ν,λ,ρ])  #:= allocates a new array
        end

        #Compare with the analytical solution
        Γ_analytical = RelativisticDynamics.christoffel(r,θ,a)
        @test isapprox(Γ,Γ_analytical)


    end

end 






@testset "Riemann symbols" begin
    

    for n in 1:1

        t = 0.0
        r = 20.0
        θ = π/4.0
        ϕ = π/6.0 
        a = 0.0

        #Calculate the metric 
        g = RelativisticDynamics.covariant_metric(r,θ,a)
        Δ = RelativisticDynamics.delta(r,a)
        g_inverse = RelativisticDynamics.contravariant_metric(g,Δ*sin(θ)^2)


        # First derivatives of the covariant metric 
        g_∂ = zeros(Float64,4,4,4) # a derivative tensor 

        g_∂[1,1,:] = jacobian(x -> RelativisticDynamics.metric_g11(x,a), [t r θ ϕ])[1]
        g_∂[2,2,:] = jacobian(x -> RelativisticDynamics.metric_g22(x,a), [t r θ ϕ])[1]
        g_∂[3,3,:] = jacobian(x -> RelativisticDynamics.metric_g33(x,a), [t r θ ϕ])[1]
        g_∂[4,4,:] = jacobian(x -> RelativisticDynamics.metric_g44(x,a), [t r θ ϕ])[1]
        g_∂[1,4,:] = jacobian(x -> RelativisticDynamics.metric_g14(x,a), [t r θ ϕ])[1]
        g_∂[4,1,:,:] = g_∂[1,4,:,:]

        # First derivatives of the contravariant metric 
        g_inverse_∂ = zeros(Float64,4,4,4) # a derivative tensor 

        g_inverse_∂[1,1,:] = jacobian(x -> RelativisticDynamics.metric_contra_g11(x,a), [t r θ ϕ])[1]
        g_inverse_∂[2,2,:] = jacobian(x -> RelativisticDynamics.metric_contra_g22(x,a), [t r θ ϕ])[1]
        g_inverse_∂[3,3,:] = jacobian(x -> RelativisticDynamics.metric_contra_g33(x,a), [t r θ ϕ])[1]
        g_inverse_∂[4,4,:] = jacobian(x -> RelativisticDynamics.metric_contra_g44(x,a), [t r θ ϕ])[1]
        g_inverse_∂[1,4,:] = jacobian(x -> RelativisticDynamics.metric_contra_g14(x,a), [t r θ ϕ])[1]
        g_inverse_∂[4,1,:,:] = g_inverse_∂[1,4,:,:]

        # Second derivative of the metric 
        g_∂∂ = zeros(Float64,4,4,4,4) # a derivative tensor
        g_∂∂[1,1,:,:] = hessian(x -> metric_g11(x,a), [t r θ ϕ])
        g_∂∂[2,2,:,:] = hessian(x -> metric_g22(x,a), [t r θ ϕ])
        g_∂∂[3,3,:,:] = hessian(x -> metric_g33(x,a), [t r θ ϕ])
        g_∂∂[4,4,:,:] = hessian(x -> metric_g44(x,a), [t r θ ϕ])
        g_∂∂[1,4,:,:] = hessian(x -> metric_g14(x,a), [t r θ ϕ])
        g_∂∂[4,1,:,:] = g_∂∂[1,4,:,:]



        @tensor begin
            Riemann[ρ,σ,μ,ν] := 0.50*(g_inverse_∂[ρ,d,μ]*(g_∂[d,ν,σ] + g_∂[d,σ,ν] - g_∂[ν,σ,d]) 
                                     +g_inverse[ρ,d]*(g_∂∂[d,ν,σ,μ]+g_∂∂[d,σ,ν,μ]-g_∂∂[ν,σ,d,μ])
                                    - g_inverse_∂[ρ,e,ν]*(g_∂[e,μ,σ] + g_∂[e,σ,μ] - g_∂[μ,σ,e]) 
                                    -g_inverse[ρ,e]*(g_∂∂[e,μ,σ,ν]+g_∂∂[e,σ,μ,ν]-g_∂∂[μ,σ,e,ν])
                                     )
        end



    end

end 



# Riemann[μ,ν,ρ,σ] := 0.50*(g_inverse_∂[μ,λ,ρ]*(g_∂[λ,ν,σ] + g_∂[λ,σ,ν] - g_∂[ν,σ,λ])
# + g_inverse[μ,λ]*(g_∂∂[λ,ν,σ,ρ] + g_∂∂[λ,σ,ν,ρ] - g_∂∂[ν,σ,λ,ρ])
# )
# -0.50*(g_inverse_∂[μ,λ,σ]*(g_∂[λ,ν,ρ] + g_∂[λ,ρ,ν] - g_∂[ν,ρ,λ])
# + g_inverse[μ,λ]*(g_∂∂[λ,ν,ρ,σ] + g_∂∂[λ,ρ,ν,σ] - g_∂∂[ν,ρ,λ,σ])
# )

# Riemann[μ,ν,ρ,σ] := 0.50*(g_inverse_∂[μ,λ,ρ]*(g_∂[λ,ν,σ] + g_∂[λ,σ,ν] - g_∂[ν,σ,λ])
                                      
# )
# -0.50*(g_inverse_∂[μ,λ,σ]*(g_∂[λ,ν,ρ] + g_∂[λ,ρ,ν] - g_∂[ν,ρ,λ])
 
#  )



#  #blob = gradient(x -> RelativisticDynamics.metric_gtt(x,θ,a), r)[1]
#         #println(blob)
   
#         gbar = zeros(Float64,4,4,4,4) # a second derivative tensor

#         t = 100.0 #arbitrary
#         ϕ = 5.0 #arbitrary 
        
            





#         @tensor begin
#             R[μ,ν,ρ,σ] := 0.50*g_inverse[μ,ρ]*(g_∂[ρ,ν,λ] +g_∂[ρ,λ,ν]-g_∂[ν,λ,ρ])  #:= allocates a new array
#         end