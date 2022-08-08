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


        r = 20.0
        θ = π/4.0 
        a = 0.0

        println(r," ",θ," ", a)

        #Calculate the metric 
        g = RelativisticDynamics.covariant_metric(r,θ,a)
        Δ = RelativisticDynamics.delta(r,a)
        g_inverse = RelativisticDynamics.contravariant_metric(g,Δ*sin(θ)^2)

        #Use automatic diff to get the gradient of the metric 
        #Only calculate r and θ derivatives
        #`gradient()` returns a tuple so [1] extracts the float 
        g_∂r = zeros(Float64,size(g))
        g_∂r[1,1] = gradient(x -> RelativisticDynamics.metric_gtt(x,θ,a), r)[1]
        g_∂r[2,2] = gradient(x -> RelativisticDynamics.metric_grr(x,θ,a), r)[1]
        g_∂r[3,3] = gradient(x -> RelativisticDynamics.metric_gθθ(x,θ,a), r)[1]
        g_∂r[4,4] = gradient(x -> RelativisticDynamics.metric_gϕϕ(x,θ,a), r)[1]
        g_∂r[4,1] = gradient(x -> RelativisticDynamics.metric_gtϕ(x,θ,a), r)[1]
        g_∂r[1,4] = g_∂r[4,1]


        g_∂θ = zeros(Float64,size(g))
        g_∂θ[1,1] = gradient(x -> RelativisticDynamics.metric_gtt(r,x,a), θ)[1]
        g_∂θ[2,2] = gradient(x -> RelativisticDynamics.metric_grr(r,x,a), θ)[1]
        g_∂θ[3,3] = gradient(x -> RelativisticDynamics.metric_gθθ(r,x,a), θ)[1]
        g_∂θ[4,4] = gradient(x -> RelativisticDynamics.metric_gϕϕ(r,x,a), θ)[1]
        g_∂θ[4,1] = gradient(x -> RelativisticDynamics.metric_gtϕ(r,x,a), θ)[1]
        g_∂θ[1,4] = g_∂θ[4,1]

        # Use these gradients to get Christoffel symbols 
        g_∂ = zeros(Float64,4,4,4)
        g_∂[:,:,2] = g_∂r
        g_∂[:,:,3] = g_∂θ
        @tensor begin
            Γ[μ,ν,λ] := 0.50*g_inverse[μ,ρ]*(g_∂[ρ,ν,λ] +g_∂[ρ,λ,ν]-g_∂[ν,λ,ρ])  #:= allocates a new array
        end

        #Compare with the analytical solution
        Γ_analytical = RelativisticDynamics.christoffel(r,θ,a)
        @test isapprox(Γ,Γ_analytical)
   
    end

end 