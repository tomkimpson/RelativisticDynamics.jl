using RelativisticDynamics
using Test


@testset "Basic run through" begin
    try
        solution,model = RelativisticDynamics.orbit()
        @test true # Completes without any errors 
    catch e
        @test  false
    end
end



@testset "Boundschecks" begin



    #Spin parameter out of bounds
    try
        solution,model = RelativisticDynamics.orbit(a=1.1)
        @test false 
    catch e
        @test  true
    end

    #Relative masses not correct
    try
        solution,model = RelativisticDynamics.orbit(mBH =1e6, mPSR=1e6)
        @test false  
    catch e
        @test  true
    end


    #Pulsar radius is unphysical
    try
        solution,model = RelativisticDynamics.orbit(rPSR=1e6)
        @test false  
    catch e
        @test  true
    end


    #Eccentricity
    try
        solution,model = RelativisticDynamics.orbit(e=1.9)
        @test false  
    catch e
        @test  true
    end

    #Orbit dir
    try
        solution,model = RelativisticDynamics.orbit(orbit_dir=1.1)
        @test false  
    catch e
        @test  true
    end


    #Incliation
    try
        solution,model = RelativisticDynamics.orbit(ι=10.0)
        @test false  
    catch e
        @test  true
    end

end