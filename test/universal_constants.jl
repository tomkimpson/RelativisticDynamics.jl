using RelativisticDynamics
using Test
using Zygote
using TensorOperations
using LinearAlgebra
using Distributions


using Combinatorics




@testset "Zero spin MPD case for mapping functions" begin
    
    NF = Float64

    m = :MPD
    a = 0.0

    for i in 1:5
        r = rand(Uniform(3.0,1e5)) 
        θ = rand(Uniform(0.0, 2.0*π))
       # P = SystemParameters(NF=NF,r=r,θ=θ,a=a,model=m) # Parameters
        
        #Mapping functions
        fr = RelativisticDynamics.mapping_f(r,a,cos(θ))
        gr = RelativisticDynamics.mapping_g(r,a)
        hr = RelativisticDynamics.mapping_h(r,a,cos(θ))
        dr = RelativisticDynamics.mapping_d(r,a,cos(θ))

        #Zero spin definitions
        @test fr == r^4
        @test gr == 0.0
        @test isapprox(hr, r*(r-2.0) + (cos(θ)^2)/(1.0 - cos(θ)^2) *(r^2-2.0*r))
        @test isapprox(dr, (r^2)*(r^2-2.0*r))

    end 
    
end


@testset "MPD Carter constant for zero inclination" begin
    
    NF = Float64

    m = :MPD
    ι= π/2.0  
    P = SystemParameters(NF=NF,ι=ι,model=m)
    
    E,L,Q = RelativisticDynamics.ELQ(P.model,P.a,P.α,P.e,P.ι,P.orbit_dir)
    @test isapprox(Q,0.0,atol=eps(NF))

  

end


@testset "MPD Angular momentum changes with radius" begin
    
    NF = Float64

    m = :MPD
    α = 1000.0
    P = SystemParameters(NF=NF,α=α,model=m)
    
    E0,L0,Q0 = RelativisticDynamics.ELQ(P.model,P.a,P.α,P.e,P.ι,P.orbit_dir)



    for i in 1:5
        α = rand(Uniform(3.0,900)) 
        P = SystemParameters(NF=NF,α=α,model=m) # Parameters
   
        E,L,Q = RelativisticDynamics.ELQ(P.model,P.a,P.α,P.e,P.ι,P.orbit_dir)
        @test L0 > L

    end 


  

end



@testset "Handles udefined models correctly" begin
    
    NF = Float64
    m = :NotMPD #some undefined model 
    P = SystemParameters(NF=NF,model=m)

    try
        C = Constants(P)
        @test false 
    catch e
        @test  true #should throw an error
    end

end



@testset "Check call of constants" begin
    
    NF = Float64
    for m in [:MPD]

        for n in 1:5

            α    = rand(Uniform(3.0,900)) 
            mBH  = rand(Uniform(1e3, 1e9))
            mPSR = rand(Uniform(1.1, 2.1))
            p0   = rand(Uniform(1e-4, 1.0))
            e    = rand(Uniform(0.01, 0.90))

            P = SystemParameters(NF=NF,α=α,mBH=mBH,mPSR=mPSR,p0=p0,e=e,model=m)

            try
                C = Constants(P)
                @test true 
            catch e
                @test false 
            end


        end 
        

    end 


end
























