using RelativisticDynamics
using Test


@testset "Basic run through" begin

        NF = Float64
        P = RelativisticDynamics.SystemParameters(NF=NF)
        C = RelativisticDynamics.Constants(P)                      
        M = RelativisticDynamics.Model(P,C)                       
        initialization = RelativisticDynamics.initial_conditions(M)
        
        try
            solution = RelativisticDynamics.timestepping(initialization, M)
            @test true # Completes without any errors 
        catch e
            @test  false
        end
end

