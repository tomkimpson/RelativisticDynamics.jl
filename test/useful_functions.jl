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
    @test K == 48.0/r^6


    fr = RelativisticDynamics.mapping_f(r,a,cos(θ))
    gr = RelativisticDynamics.mapping_g(r,a)
    hr = RelativisticDynamics.mapping_h(r,a,cos(θ))
    dr = RelativisticDynamics.mapping_d(r,a,cos(θ))

    @test fr == r^4
    @test gr == 0.0
    @test hr == r*(r-2.0) + (cos(θ)^2)/(1.0 - cos(θ)^2) *(r^2-2.0*r)
    @test dr == (r^2)*(r^2-2.0*r)


    

end












