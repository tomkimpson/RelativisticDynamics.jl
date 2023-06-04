@testset "Raise and lower vector indices" begin
    
    for n in 1:5
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


@testset "Levi civita is antisymmetric" begin
    
    for n in 1:5
        r = rand(Uniform(3.0,1e5))      # Radial coordinate. 3.0 as rough lower limit of an event horizon
        θ = rand(Uniform(0.0, 2.0*π))
        a = rand(Uniform(0.0, 1.0))


        xvector=[0.0,r,θ,0.0] 
        g  = RelativisticDynamics.covariant_metric(xvector,a)

        levi = RelativisticDynamics.levi_civita_tensor(g)  #This is the fully contravariant Levi Civita tensor 
        @test isapprox(permutedims(levi, (2,1,3,4)),-levi) #exchange of just 1 and 2

    end 
    
end













