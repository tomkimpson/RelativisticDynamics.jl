@testset "Kerr functions for a=0" begin
    
    r = rand(Uniform(3.0,1e5))      # Radial coordinate. 3.0 as rough lower limit of an event horizon
    θ = rand(Uniform(0.0, 2.0*π))
    a = 0.0


    Δ = RelativisticDynamics.delta(r,a)
    Σ =  RelativisticDynamics.sigma(r,θ,a)
 
    #This is what they should be if a=0.0
    @test Δ == r^2 - 2.0*r
    @test Σ == r^2 

    K = RelativisticDynamics.Kretschmann_scalar(r,θ,a)
    @test isapprox(K,48.0/r^6,atol=eps(Float64))

end


@testset "Raise and lower vector indices" begin
    
    for n in 1:1
        r = rand(Uniform(3.0,1e5))      # Radial coordinate. 3.0 as rough lower limit of an event horizon
        θ = rand(Uniform(0.0, 2.0*π))
        a = rand(Uniform(0.0, 1.0))


        xvector=[0.0,r,θ,0.0] 
        g  = RelativisticDynamics.covariant_metric(xvector,a)

        contravariant_vector = [1.0,2.0,3.0,4.0]
        covariant_vector = RelativisticDynamics.convert_to_covariant(g,contravariant_vector)

        ginv  = RelativisticDynamics.contravariant_metric(xvector,a) #convert_to_covariant is poorly named. Really raises/lowers depending on input metric
        new_contravariant_vector = RelativisticDynamics.convert_to_covariant(ginv,covariant_vector)

        @test isapprox(contravariant_vector,new_contravariant_vector)
    end 
    
end











