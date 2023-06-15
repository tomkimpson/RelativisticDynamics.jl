@testset "MPD Carter constant for zero inclination" begin
    
    NF = Float64

    
    ι= π/2.0  
    P = RelativisticDynamics.SystemParameters(NF=NF,ι=ι)
    
    E,L,Q = RelativisticDynamics.ELQ(P.a,P.α,P.e,P.ι,P.orbit_dir)
    @test isapprox(Q,0.0,atol=eps(NF))

  

end


@testset "MPD Angular momentum changes with radius" begin
    
    NF = Float64

    α = 1000.0
    P = RelativisticDynamics.SystemParameters(NF=NF,α=α)
    
    E0,L0,Q0 = RelativisticDynamics.ELQ(P.a,P.α,P.e,P.ι,P.orbit_dir)

    for i in 1:5
        α = rand(Uniform(3.0,900)) 
        P = RelativisticDynamics.SystemParameters(NF=NF,α=α) # Parameters
   
        E,L,Q = RelativisticDynamics.ELQ(P.a,P.α,P.e,P.ι,P.orbit_dir)
        @test L0 > L

    end 

end

@testset "Check basic call of constants" begin
    
    NF = Float64

    for n in 1:5

        α    = rand(Uniform(3.0,900)) 
        mBH  = rand(Uniform(1e3, 1e9))
        mPSR = rand(Uniform(1.1, 2.1))
        p0   = rand(Uniform(1e-4, 1.0))
        e    = rand(Uniform(0.01, 0.90))

        P = RelativisticDynamics.SystemParameters(NF=NF,α=α,mBH=mBH,mPSR=mPSR,p0=p0,e=e)
        
        try
            C = RelativisticDynamics.Constants(P)
            @test true 
        catch e
            @test false 
        end


    end 
    

end


