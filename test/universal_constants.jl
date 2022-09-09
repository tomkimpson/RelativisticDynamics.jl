using RelativisticDynamics
using Test
using Zygote
using TensorOperations
using LinearAlgebra
using Distributions


using Combinatorics

@testset "No spherical photon orbits in Schwarzchild" begin
    
  
    NF = Float64

    m = :SphericalPhoton
    r = rand(Uniform(3.0,1e5))      
    θ = rand(Uniform(0.0, 2.0*π))
    a = 0.0

    P = SystemParameters(NF=NF,r=r,θ=θ,a=a,model=m) # Parameters
   
   #a=0 should throw error for spherical photon orbits
    try
        L,Q = RelativisticDynamics.LQ(P)
    catch e
        @test true 
    end

end



@testset "No θ dependence for spherical photon orbits L/Q" begin
    
    NF = Float64

    m = :SphericalPhoton
    r = rand(Uniform(3.0,1e5)) 
    a = rand(Uniform(-0.99, 0.99))
    θ = rand(Uniform(0.0, 2.0*π))
    P = SystemParameters(NF=NF,r=r,θ=θ,a=a,model=m) # Parameters

    L0,Q0 = RelativisticDynamics.LQ(P)

    for i in 1:5

        θ = rand(Uniform(0.0, 2.0*π))
        P = SystemParameters(NF=NF,r=r,θ=θ,a=a,model=m) # Parameters
        L,Q = RelativisticDynamics.LQ(P)

        @test L == L0 && Q == Q0

    end 

end


@testset "Zero spin MPD case for mapping functions" begin
    
    NF = Float64

    m = :MPD
    a = 0.0

    for i in 1:5
        r = rand(Uniform(3.0,1e5)) 
        θ = rand(Uniform(0.0, 2.0*π))
        P = SystemParameters(NF=NF,r=r,θ=θ,a=a,model=m) # Parameters
        
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
    
    L,Q = RelativisticDynamics.LQ(P)
    @test isapprox(Q,0.0,atol=eps(NF))

  

end


@testset "MPD Angular momentum changes with radius" begin
    
    NF = Float64

    m = :MPD
    α = 1000.0
    P = SystemParameters(NF=NF,α=α,model=m)
    
    L0,Q0 = RelativisticDynamics.LQ(P)



    for i in 1:5
        α = rand(Uniform(3.0,900)) 
        θ = rand(Uniform(0.0, 2.0*π))
        P = SystemParameters(NF=NF,α=α,θ=θ,model=m) # Parameters
   
        L,Q = RelativisticDynamics.LQ(P)
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

    m = :MPD

    for m in [:SphericalPhoton, :MPD]

        for n in 1:5

            α    = rand(Uniform(3.0,900)) 
            θ    = rand(Uniform(0.0, 2.0*π))
            mBH  = rand(Uniform(1e3, 1e9))
            mPSR = rand(Uniform(1.1, 2.1))
            p0   = rand(Uniform(1e-4, 1.0))
            e    = rand(Uniform(0.01, 0.90))

            P = SystemParameters(NF=NF,α=α,θ=θ,mBH=mBH,mPSR=mPSR,p0=p0,e=e,model=m)

            try
                C = Constants(P)
                @test true 
            catch e
                @test false 
            end


        end 
        

    end 


end
























### TBD - confirm that covariant eps = -contravariant eps 



# r = rand(Uniform(3.0,1e5))      # Radial coordinate. 3.0 as rough lower limit of an event horizon
# θ = rand(Uniform(0.0, 2.0*π))
# a = 0.0

# Σ = RelativisticDynamics.sigma(r,θ,a)
# g = RelativisticDynamics.covariant_minkowski() ##RelativisticDynamics.covariant_metric([rand(Float64, 1),r,θ,rand(Float64, 1)],a)



# #metric_trace =-sin(θ)^2*Σ^2
# metric_trace = 2.0
# permutation_tensor = C.ϵ/sqrt(abs(metric_trace)) 

# # Turn the contravariant permutation_tensor to a fully covariant form
# @tensor begin
#     eps[ρ,σ,μ,ν] := g[ρ,α]*g[σ,β]*g[μ,λ]*g[ν,γ]*permutation_tensor[α,β,λ,γ] 
# end

# for i in 1:4
#     for j in 1:4
#         for k in 1:4
#             for l in 1:4
#                 println(i,j,k,l,"   ",permutation_tensor[i,j,k,l], "  ",eps[i,j,k,l] )
#             end
#         end 
#     end
# end 

# permutation_vector = [1,3,2,4]
# println(levicivita(permutation_vector))










# @testset "Constants of motion" begin
    
  
#     NF = Float64

#     for m in [:SphericalPhoton]

#         P = SystemParameters(NF=NF,model=m) # Parameters

#         #Create some random, Schwarzchild coordiantes
#         r = rand(Uniform(3.0,1e5))      # Radial coordinate. 3.0 as rough lower limit of an event horizon
#         θ = rand(Uniform(0.0, 2.0*π))
#         a = 0.0




#     end 

#     #C = Constants(P)                                   # Constants


   

#     fr = RelativisticDynamics.mapping_f(r,a,cos(θ))
#     gr = RelativisticDynamics.mapping_g(r,a)
#     hr = RelativisticDynamics.mapping_h(r,a,cos(θ))
#     dr = RelativisticDynamics.mapping_d(r,a,cos(θ))

#     @test fr == r^4
#     @test gr == 0.0
#     @test hr == r*(r-2.0) + (cos(θ)^2)/(1.0 - cos(θ)^2) *(r^2-2.0*r)
#     @test dr == (r^2)*(r^2-2.0*r)



#     L,Q = RelativisticDynamics.LQ(P)

#     println(L)
#     println(Q)


# end









