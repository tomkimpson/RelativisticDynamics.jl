@testset "Kerr functions" begin
    

    r = 10.0
    θ = π/6.0   
    a = 0.0


    Δ = RelativisticDynamics.delta(r,a)
    Σ =  RelativisticDynamics.sigma(r,θ,a)

    @test Δ == r^2 - 2.0*r
    @test Σ == r^2 





end












