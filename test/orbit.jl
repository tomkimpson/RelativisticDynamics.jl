@testset "Basic run through" begin
    try
        solution,model = RelativisticDynamics.orbit()
        @test true # The default completes without any errors 
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
        @test e == ErrorException("Spin parameter a is out of bounds")
    end

    #Relative masses not correct
    try
        solution,model = RelativisticDynamics.orbit(mBH =1e6, mPSR=1e6)
        @test false  
    catch e
        @test e == ErrorException("Mass ratio is too small")
    end


    #Pulsar radius is unphysical
    try
        solution,model = RelativisticDynamics.orbit(rPSR=1e6)
        @test false  
    catch e
        @test e == ErrorException("Pulsar radius is unphysical")
    end


    #Eccentricity
    try
        solution,model = RelativisticDynamics.orbit(e=1.9)
        @test false  
    catch e
        @test e == ErrorException("Eccentricity is outside range")
    end

    #Orbit dir must be 1/-1
    try
        solution,model = RelativisticDynamics.orbit(orbit_dir=2)
        @test false  
    catch e
        @test e == ErrorException("Orbit direction must be plus or minus 1")
    end

    #Orbit dir must be integer
    try
        solution,model = RelativisticDynamics.orbit(orbit_dir=1.1)
        @test false  
    catch e
        @test e == InexactError(:Int64, Int64, 1.1)
    end 

    #Incliation
    try
        solution,model = RelativisticDynamics.orbit(ι=10.0)
        @test false  
    catch e
        @test e == ErrorException("ι is outside of allowed range")
    end

end