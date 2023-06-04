
@testset "Check all number formats are consistent" begin

    for NF in [Float64,Float32]
        solution,model = RelativisticDynamics.orbit(NF=NF)
        @test eltype(solution) == NF

    end

end 