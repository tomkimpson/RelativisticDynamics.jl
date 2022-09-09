
using RelativisticDynamics
using Test
using Zygote
using TensorOperations
using LinearAlgebra
using Distributions



@testset "Each model maps to correct initialization function" begin
    
    NF = Float64

    m =:SphericalPhoton
    P = SystemParameters(NF=NF,model=m)
    C = Constants(P)
    M = RelativisticDynamics.Model(P,C) 
    prognostic_vars = RelativisticDynamics.initial_conditions(M)

    @test prognostic_vars.svector == [0,0,0,0] # For spherical photon orbits we do not initialise spin initial conditions 
    @test prognostic_vars.pvector[1] == 0.0 
    @test prognostic_vars.pvector[2] == 0.0
    @test prognostic_vars.pvector[3] != 0.0
    @test prognostic_vars.pvector[4] == 0.0




    m =:MPD
    P = SystemParameters(NF=NF,model=m)
    C = Constants(P)
    M = RelativisticDynamics.Model(P,C) 
    prognostic_vars = RelativisticDynamics.initial_conditions(M)

    @test prognostic_vars.svector != [0,0,0,0] # For spherical photon orbits we do not initialise spin initial conditions 
    @test prognostic_vars.pvector[1] != [0,0,0,0]

    # 1. Four- position
    @test prognostic_vars.xvector == [0.0,P.α,P.θ,P.ϕ] 

    Stensor = RelativisticDynamics.spintensor(prognostic_vars.xvector,prognostic_vars.pvector,prognostic_vars.svector,P.a,C.m0,C.ϵ)


    @tensor begin
        blob[μ] := Stensor[μ,ν]*prognostic_vars.pvector[ν]
    end

    display(blob)

end











