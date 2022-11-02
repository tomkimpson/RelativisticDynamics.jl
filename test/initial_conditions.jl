
using RelativisticDynamics
using Test
using Zygote
using TensorOperations
using LinearAlgebra
using Distributions


@testset "Initial conditions runs OK for default initialization" begin
    
    NF = Float64
    P = SystemParameters(NF=NF)
    C = Constants(P)                      
    M = RelativisticDynamics.Model(P,C)                       

    try
        initialization = RelativisticDynamics.initial_conditions(M)
        @test true 
    catch e
        @test  false 
    end
    

end


@testset "Initial values of xvector are always as expected" begin
    
    NF = Float64


    for n in 1:5
        α    = rand(Uniform(10.0,500)) 
        e    = rand(Uniform(0.01, 0.90))

        P = SystemParameters(NF=NF,α=α,e=e)
        C = Constants(P)                                   
        M = RelativisticDynamics.Model(P,C)                       


        initialization = RelativisticDynamics.initial_conditions(M)

        @test initialization.xvector[1] == 0.0
        @test initialization.xvector[2] == α
        @test initialization.xvector[3] == π/2.0
        @test initialization.xvector[4] == 0.0
    end

end



@testset "4-velocities are reasonable / not complex" begin
    
    NF = Float64


    for n in 1:5
        α    = rand(Uniform(3.0,900)) 
        ι    = rand(Uniform(0.0,π/2.0)) 
        e    = rand(Uniform(0.01, 0.90))

        P = SystemParameters(NF=NF,α=α,e=e,ι=ι)
        C = Constants(P)                                   
        M = RelativisticDynamics.Model(P,C)                       


        initialization = RelativisticDynamics.initial_conditions(M)

        @test typeof(initialization.pvector[1])<:Real
        @test typeof(initialization.pvector[2])<:Real
        @test typeof(initialization.pvector[3])<:Real
        @test typeof(initialization.pvector[4])<:Real
    end

end


@testset "Orbit in the plane" begin
    
    NF = Float64


    for n in 1:5
        α    = rand(Uniform(3.0,900)) 
        ι    = π/2.0
        e    = rand(Uniform(0.01, 0.90))

        P = SystemParameters(NF=NF,α=α,e=e,ι=ι)
        C = Constants(P)                                   
        M = RelativisticDynamics.Model(P,C)                       


        initialization = RelativisticDynamics.initial_conditions(M)
        @test initialization.pvector[3] == 0.0 #no theta dot
    end

end



@testset "TD spin condition" begin
    
    NF = Float64


    for n in 1:5
        α    = rand(Uniform(3.0,900)) 
        ι    = rand(Uniform(0.0,π/2.0))
        e    = rand(Uniform(0.01, 0.90))
        a    = rand(Uniform(0.01, 0.99))


        P = SystemParameters(NF=NF,α=α,e=e,ι=ι,a=a)
        C = Constants(P)                                   
        M = RelativisticDynamics.Model(P,C)                       


        initialization = RelativisticDynamics.initial_conditions(M)


        g  = RelativisticDynamics.covariant_metric(initialization.xvector,a)


        levi = RelativisticDynamics.permutation_tensor(g,C.ϵ)  #This is the fully contravariant Levi Civita tensor 
        spin_tensor = RelativisticDynamics.spintensor(levi,initialization.pvector,initialization.svector,C.m0) #the fully contravariant spin tensor s^{ab}

        pvector_covar = RelativisticDynamics.convert_to_covariant(g,initialization.pvector)

        @tensor begin
            val[μ] := spin_tensor[μ,ν]*pvector_covar[ν]
       end
        
       @test isapprox(val,[0.0,0.0,0.0,0.0],atol=eps(NF))

    end

end









