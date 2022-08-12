@testset "Trace of the metric" begin
    
    for n in 1:3

        #Get some coordiantes at random 
        r = rand(Uniform(3.0,1e5))      # Radial coordinate. 3.0 as rough lower limit of an event horizon
        θ = rand(Uniform(0.0, 2.0*π))   # Polar coord
        a = rand(Uniform(-0.99, 0.99))  # Spin parameter
        
        #Calculate the metric 
        coords = [0.0,r,θ,0.0]
        g = RelativisticDynamics.covariant_metric(coords,a)
        g_inverse = RelativisticDynamics.contravariant_metric(coords,a)
        #Check metric trace is OK 
        @tensor begin
        δ[a,c] := g[a,b] * g_inverse[b,c]  #:= allocates a new array
        end
        @test tr(δ)==4.0 
    end

end




@testset "Christoffel tensor components" begin
    

    for n in 1:3

        #Get some coordiantes at random 
        t = rand(Uniform(3.0,1e5))       # Time coordinate - arbitratry since metric is time-independent
        r = rand(Uniform(3.0,1e5))       # Radial coordinate. 3.0 as rough lower limit of an event horizon
        θ = rand(Uniform(0.0, 2.0*π))    # Polar coordinate
        ϕ = rand(Uniform(0.0, 2.0*π))    # Azimuth coordinate - arbitratry since metric is axisymmetric
        a = rand(Uniform(-0.99, 0.99))   # Spin parameter

        coords = [t,r,θ,ϕ]
        
        #Calculate the metric 
        g = RelativisticDynamics.covariant_metric(coords,a)
        g_inverse = RelativisticDynamics.contravariant_metric(coords,a)


        #Use automatic diff to get the first derivatives of the metric  of the metric 
        #Only calculate r and θ derivatives
        #`jacobian()` returns a tuple so [1] extracts the vector of floats 
        g_∂ = zeros(Float64,4,4,4) # a derivative tensor 

        g_∂[1,1,:] = jacobian(x -> RelativisticDynamics.metric_g11(x,a), [t r θ ϕ])[1]
        g_∂[2,2,:] = jacobian(x -> RelativisticDynamics.metric_g22(x,a), [t r θ ϕ])[1]
        g_∂[3,3,:] = jacobian(x -> RelativisticDynamics.metric_g33(x,a), [t r θ ϕ])[1]
        g_∂[4,4,:] = jacobian(x -> RelativisticDynamics.metric_g44(x,a), [t r θ ϕ])[1]
        g_∂[1,4,:] = jacobian(x -> RelativisticDynamics.metric_g14(x,a), [t r θ ϕ])[1]

        @tensor begin
            Γ[μ,ν,λ] := 0.50*g_inverse[μ,ρ]*(g_∂[ρ,ν,λ] +g_∂[ρ,λ,ν]-g_∂[ν,λ,ρ])  #:= allocates a new array
        end

        #Compare with the analytical solution
        Γ_analytical = RelativisticDynamics.christoffel(r,θ,a)
        @test isapprox(Γ,Γ_analytical)


    end

end 


@testset "Riemann tensor components" begin
    

    for n in 1:5

        #Get some coordiantes at random 
        t = rand(Uniform(3.0,1e5))       # Time coordinate - arbitratry since metric is time-independent
        r = rand(Uniform(3.0,1e3))       # Radial coordinate. 3.0 as rough lower limit of an event horizon
        θ = rand(Uniform(0.0, 2.0*π))    # Polar coordinate
        ϕ = rand(Uniform(0.0, 2.0*π))    # Azimuth coordinate - arbitratry since metric is axisymmetric
        a = rand(Uniform(-0.99, 0.99))   # Spin parameter

 
        coords = [t,r,θ,ϕ]

        #Calculate the metric 
        g = RelativisticDynamics.covariant_metric(coords,a)
        g_inverse = RelativisticDynamics.contravariant_metric(coords,a)

        # First derivatives of the covariant metric 
        g_∂ = zeros(Float64,4,4,4) # a derivative tensor 

        g_∂[1,1,:] = jacobian(x -> RelativisticDynamics.metric_g11(x,a), [t r θ ϕ])[1]
        g_∂[2,2,:] = jacobian(x -> RelativisticDynamics.metric_g22(x,a), [t r θ ϕ])[1]
        g_∂[3,3,:] = jacobian(x -> RelativisticDynamics.metric_g33(x,a), [t r θ ϕ])[1]
        g_∂[4,4,:] = jacobian(x -> RelativisticDynamics.metric_g44(x,a), [t r θ ϕ])[1]
        g_∂[1,4,:] = jacobian(x -> RelativisticDynamics.metric_g14(x,a), [t r θ ϕ])[1]
        g_∂[4,1,:,:] = g_∂[1,4,:,:]

        # First derivatives of the contravariant metric 
        g_inverse_∂ = zeros(Float64,4,4,4) # a derivative tensor 

        g_inverse_∂[1,1,:] = jacobian(x -> RelativisticDynamics.metric_contra_g11(x,a), [t r θ ϕ])[1]
        g_inverse_∂[2,2,:] = jacobian(x -> RelativisticDynamics.metric_contra_g22(x,a), [t r θ ϕ])[1]
        g_inverse_∂[3,3,:] = jacobian(x -> RelativisticDynamics.metric_contra_g33(x,a), [t r θ ϕ])[1]
        g_inverse_∂[4,4,:] = jacobian(x -> RelativisticDynamics.metric_contra_g44(x,a), [t r θ ϕ])[1]
        g_inverse_∂[1,4,:] = jacobian(x -> RelativisticDynamics.metric_contra_g14(x,a), [t r θ ϕ])[1]
        g_inverse_∂[4,1,:,:] = g_inverse_∂[1,4,:,:]

        # Second derivative of the metric 
        g_∂∂ = zeros(Float64,4,4,4,4) # a derivative tensor
        g_∂∂[1,1,:,:] = hessian(x -> RelativisticDynamics.metric_g11(x,a), [t r θ ϕ])
        g_∂∂[2,2,:,:] = hessian(x -> RelativisticDynamics.metric_g22(x,a), [t r θ ϕ])
        g_∂∂[3,3,:,:] = hessian(x -> RelativisticDynamics.metric_g33(x,a), [t r θ ϕ])
        g_∂∂[4,4,:,:] = hessian(x -> RelativisticDynamics.metric_g44(x,a), [t r θ ϕ])
        g_∂∂[1,4,:,:] = hessian(x -> RelativisticDynamics.metric_g14(x,a), [t r θ ϕ])
        g_∂∂[4,1,:,:] = g_∂∂[1,4,:,:]


        # e.g. https://math.stackexchange.com/questions/628863/riemann-tensor-in-terms-of-the-metric-tensor
        #This is R^{\rho}_{\sigma \mu \nu} 
        @tensor begin
            Riemann[ρ,σ,μ,ν] := 
                                0.50*(g_inverse_∂[ρ,d,μ]*(g_∂[d,ν,σ] + g_∂[d,σ,ν] - g_∂[ν,σ,d])  + 
                                      g_inverse[ρ,d]*(g_∂∂[d,ν,σ,μ]+g_∂∂[d,σ,ν,μ]-g_∂∂[ν,σ,d,μ]) -
                                      g_inverse_∂[ρ,e,ν]*(g_∂[e,μ,σ] + g_∂[e,σ,μ] - g_∂[μ,σ,e])  -
                                      g_inverse[ρ,e]*(g_∂∂[e,μ,σ,ν]+g_∂∂[e,σ,μ,ν]-g_∂∂[μ,σ,e,ν])
                                     ) +
                                0.25*g_inverse[ρ,f]*(g_∂[f,μ,λ] + g_∂[f,λ,μ] - g_∂[μ,λ,f])*g_inverse[λ,h]*(g_∂[h,ν,σ] + g_∂[h,σ,ν] - g_∂[ν,σ,h]) -
                                0.25*g_inverse[ρ,i]*(g_∂[i,ν,λ] + g_∂[i,λ,ν] - g_∂[ν,λ,i])*g_inverse[λ,k]*(g_∂[k,μ,σ] + g_∂[k,σ,μ] - g_∂[μ,σ,k])
                                
        end



        #Compare with the analytical solution
        Riemann_analytical = RelativisticDynamics.riemann(r,θ,a)
        @test isapprox(Riemann,Riemann_analytical)


        #Also double check the covariant form
        @tensor begin
            Riemann_covar[μ,ν,ρ,σ] := g[μ,λ]*Riemann[λ,ν,ρ,σ]
        end
        #The purely covariant version of the Riemann tensor 
        @tensor begin
            Riemann_covar_analytical[μ,ν,ρ,σ] := g[μ,λ]*Riemann_analytical[λ,ν,ρ,σ]
        end
        @test isapprox(Riemann_covar,Riemann_covar_analytical)

    end

end 


@testset "Riemann tensor reduces to schwarzchild solution" begin

    for n in 1:5

        #Get some coordiantes at random 
        t = rand(Uniform(3.0,1e5))       # Time coordinate - arbitratry since metric is time-independent
        r = rand(Uniform(3.0,1e3))       # Radial coordinate. 3.0 as rough lower limit of an event horizon
        θ = rand(Uniform(0.0, 2.0*π))    # Polar coordinate
        ϕ = rand(Uniform(0.0, 2.0*π))    # Azimuth coordinate - arbitratry since metric is axisymmetric
        a = 0.0                          # Schwarzchild

 
        coords = [t,r,θ,ϕ]

        #Calculate the metric 
        g = RelativisticDynamics.covariant_metric(coords,a)
        g_inverse = RelativisticDynamics.contravariant_metric(coords,a)



        #Get the Riemann tensor solution, analytically
        Riemann_analytical = RelativisticDynamics.riemann(r,θ,a)
        

        #Get the schwazchild solution
        Riemann_covar_schwarzchild = RelativisticDynamics.schwarzchild_covariant_riemann(r,θ,a)

        #The contra/covar version of the Riemann tensor for schwarzchild
        @tensor begin
            Riemann_schwarzchild[μ,ν,ρ,σ] := g_inverse[μ,λ]*Riemann_covar_schwarzchild[λ,ν,ρ,σ]
        end

        #Compare 
        @test isapprox(Riemann_analytical,Riemann_schwarzchild)


    end

end 













