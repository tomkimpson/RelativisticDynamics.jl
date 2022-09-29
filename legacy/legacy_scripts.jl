# @testset "No spherical photon orbits in Schwarzchild" begin
    
  
#     NF = Float64

#     m = :SphericalPhoton
#     r = rand(Uniform(3.0,1e5))      
#     θ = rand(Uniform(0.0, 2.0*π))
#     a = 0.0

#     P = SystemParameters(NF=NF,r=r,θ=θ,a=a,model=m) # Parameters
   
#    #a=0 should throw error for spherical photon orbits
#     try
#         L,Q = RelativisticDynamics.LQ(P)
#     catch e
#         @test true 
#     end

# end



# @testset "No θ dependence for spherical photon orbits L/Q" begin
    
#     NF = Float64

#     m = :SphericalPhoton
#     r = rand(Uniform(3.0,1e5)) 
#     a = rand(Uniform(-0.99, 0.99))
#     θ = rand(Uniform(0.0, 2.0*π))
#     P = SystemParameters(NF=NF,r=r,θ=θ,a=a,model=m) # Parameters

#     L0,Q0 = RelativisticDynamics.LQ(P)

#     for i in 1:5

#         θ = rand(Uniform(0.0, 2.0*π))
#         P = SystemParameters(NF=NF,r=r,θ=θ,a=a,model=m) # Parameters
#         L,Q = RelativisticDynamics.LQ(P)

#         @test L == L0 && Q == Q0

#     end 

# end















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





# function LQ(P::SystemParameters)


#     @unpack r,a,α,e,ι,orbit_dir = P

#     if P.model == :SphericalPhoton 

#         if a == 0.0
#             println("Spherical Photon orbits do not exist for a=0. Try a non-zero value")
#             return
#         else
#             E = 1.0
#             L = - (r^3 - 3.0*r^2 + a^2*r + a^2) / (a*(r - 1.0))
#             Q = -(r^3 * (r^3 - 6.0*r^2 + 9.0 * r - 4.0*a^2))/(a^2 * (r - 1.0)^2)
#         end

#     elseif P.model == :RayTracing 


#         E = 0.0
#         L = 0.0
#         Q = 0.0

        

#         E = 0.0
#         L = 0.0
#         Q = 0.0

#     elseif P.model == :MPD


#         #ELQ circular 

#         N = 1-3/r






#         # #Some needed quantities
#         # r_periapsis= α*(1-e)
#         # println("PERIAPSIS")
#         # println(r_periapsis)

#         # zminus = cos(ι)
#         # f1 = mapping_f(r_periapsis,a,zminus)
#         # g1 = mapping_g(r_periapsis,a)
#         # h1 = mapping_h(r_periapsis,a,zminus)
#         # d1 = mapping_d(r_periapsis,a,zminus)
    
   
#         # #Angular momentum and Carter constant, taking E=1
#         # L = -g1/h1 + orbit_dir*sqrt(g1^2 + (f1-d1)*h1)/h1
#         # Q = zminus^2 * (L^2 /(1.0 - zminus^2))

    

#     else
#         println(P.model)
#         println("That model selection is not defined")
#         println("Please choose one of: SphericalPhoton, MPD ")
#         return
#     end 

#     return E,L,Q


# end 










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
        TD[μ] := Stensor[μ,ν]*prognostic_vars.pvector[ν] #Tulczyjew Dixon condition 
    end
    for i in 1:4
        @test isapprox(TD[i],0.0,atol=eps(NF))
    end 


end






function spherical_photon_hamiltonian!(du,u,p,λ)

    t,r,θ,ϕ,pᵗ,pʳ,pᶿ,pᵠ = u
    L,a = p


    # Define useful functions 
    Σ = sigma(r,θ,a)
    Δ = delta(r,a)

    #Position 
    du[1] = 0.0
    du[2] = 0.0
    du[3] = pᶿ/Σ
    du[4] = (2.0*r*a + (Σ - 2.0*r)*L/sin(θ)^2) / (Σ*Δ)

    #Momentum 
    du[5] = 0.0
    du[6] = 0.0
    du[7] = sin(θ)*cos(θ)*(L^2/sin(θ)^4 -a^2)/Σ
    du[8] = 0.0

    nothing #function returns nothing


end 





# # Define the ODE to be solved
# if M.parameters.model == :SphericalPhoton      
#     params = [L,a]    
#     u = vcat(X.xvector,X.pvector)
#     f = spherical_photon_hamiltonian!     
# elseif M.parameters.model == :MPD
#     params = [a,m0,ϵ]
#     u = vcat(X.xvector,X.pvector,X.svector)
#     #f = MPD!
#     f = geodesic!
# end 





## Spherical Photon Orbits

Spherical photon orbits are also included, mainly as a test case/toy model to ensure the code machinery works as expected.

The general description of the equations of motion is described in [Teo 2003](https://ui.adsabs.harvard.edu/abs/2003GReGr..35.1909T/abstract) which have been reframed here via a Hamiltonian formalism c.f. [Pu et al. 2016](https://arxiv.org/abs/1601.02063)

