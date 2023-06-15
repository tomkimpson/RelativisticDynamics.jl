using Tullio
using Distributions
using LinearAlgebra
using Zygote
@testset "Trace of the metric" begin
    
    for n in 1:5

        #Get some coordiantes at random 
        r = rand(Uniform(3.0,1e5))      # Radial coordinate. 3.0 as rough lower limit of an event horizon
        θ = rand(Uniform(0.0, 2.0*π))   # Polar coord
        a = rand(Uniform(-0.99, 0.99))  # Spin parameter
        
        #Calculate the metric 
        coords = [0.0,r,θ,0.0]
        g = RelativisticDynamics.covariant_metric(coords,a)
        g_inverse = RelativisticDynamics.contravariant_metric(coords,a)
        
        #Check metric trace is OK 
        @tullio δ[a,c] := g[a,b] * g_inverse[b,c]  #:= allocates a new array
        @test isapprox(tr(δ),4.0)
    end

end



@testset "Christoffel tensor components" begin
    

    for n in 1:5

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

        #and the metric Jacobian using Zygote
        g_jacobian = jacobian(x -> RelativisticDynamics.covariant_metric(x,a), [t r θ ϕ])[1]
        g_jacobian = reshape(g_jacobian,(4,4,4)) # reshape

        #Christoffel
        @tullio Γ[μ,ν,λ] := 0.50*g_inverse[μ,ρ]*(g_jacobian[ρ,ν,λ] +g_jacobian[ρ,λ,ν]-g_jacobian[ν,λ,ρ])  #:= allocates a new array

        #Compare with the analytical solution
        Γ_analytical = RelativisticDynamics.christoffel(coords,a)
        @test isapprox(Γ,Γ_analytical)
    end
end 


@testset "Riemann tensor reduces to schwarzchild solution" begin



    """
    R = schwarzchild_covariant_riemann(coords,a)
    Special case - the fully covariant components of the Riemann tensor for schwarzchild metric
    Used for testing 
    """
    function schwarzchild_covariant_riemann(coords,a)

        r = coords[2]
        θ = coords[3]
    
        Rtensor = zeros(typeof(a),4,4,4,4)

        Rtensor[1,2,1,2] = -2.0/r^3
        Rtensor[1,3,1,3] = (r-2.0)/r^2
        Rtensor[1,4,1,4] = (r-2.0)*sin(θ)^2/r^2
        Rtensor[2,3,2,3] = -1.0/(r-2.0)
        Rtensor[2,4,2,4] = -sin(θ)^2/(r-2.0)
        Rtensor[3,4,3,4] = r*2.0*sin(θ)^2


        #Symmetries
        Rtensor[1,2,2,1] = -Rtensor[1,2,1,2]
        Rtensor[1,3,3,1] = -Rtensor[1,3,1,3]
        Rtensor[1,4,4,1] = -Rtensor[1,4,1,4]
        Rtensor[2,3,3,2]=  -Rtensor[2,3,2,3]
        Rtensor[2,4,4,2] = -Rtensor[2,4,2,4]
        Rtensor[3,4,4,3] = -Rtensor[3,4,3,4]

        Rtensor[2,1,1,2] = -Rtensor[1,2,1,2]
        Rtensor[3,1,1,3] = -Rtensor[1,3,1,3]
        Rtensor[4,1,1,4] = -Rtensor[1,4,1,4]
        Rtensor[3,2,2,3]=  -Rtensor[2,3,2,3]
        Rtensor[4,2,2,4] = -Rtensor[2,4,2,4]
        Rtensor[4,3,3,4] = -Rtensor[3,4,3,4]

        Rtensor[2,1,2,1] = Rtensor[1,2,1,2]
        Rtensor[3,1,3,1] = Rtensor[1,3,1,3] 
        Rtensor[4,1,4,1] = Rtensor[1,4,1,4] 
        Rtensor[3,2,3,2] = Rtensor[2,3,2,3] 
        Rtensor[4,2,4,2] = Rtensor[2,4,2,4] 
        Rtensor[4,3,4,3] = Rtensor[3,4,3,4] 

        return Rtensor



    end

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
        Riemann_analytical = RelativisticDynamics.riemann(coords,a)
        

        #Get the schwazchild solution
        Riemann_covar_schwarzchild = schwarzchild_covariant_riemann(coords,a)

        #The contra/covar version of the Riemann tensor for schwarzchild
        @tullio Riemann_schwarzchild[μ,ν,ρ,σ] := g_inverse[μ,λ]*Riemann_covar_schwarzchild[λ,ν,ρ,σ]
       # end

        #Compare 
        @test isapprox(Riemann_analytical,Riemann_schwarzchild)

       
    end

end 

@testset "Kretschman scalar" begin
    

    """
    K =  Kretschmann_scalar(r,θ,a)
    Kretschman scalar for the Kerr metric
    """
    function Kretschmann_scalar(r,θ,a)
        Σ = r^2 + a^2 * cos(θ)^2
        return 48.0*(2.0*r^2-Σ)*(Σ^2-16.0*r^2*a^2*cos(θ)^2)/Σ^6
    end 


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

       
        #Rieman tensor, mixed indices
        Riemann = RelativisticDynamics.riemann(coords,a)
      

        #Fully covariant form 
        @tullio Riemann_covar[μ,ν,ρ,σ] := g[μ,λ]*Riemann[λ,ν,ρ,σ]
        
        #Fully contravariant form 
        @tullio Riemann_contra[μ,ν,α,β] := g_inverse[μ,i]*g_inverse[ν,j]*g_inverse[α,k]*g_inverse[β,l]*Riemann_covar[i,j,k,l]
        

        #Kretschman scalar 
        @tullio Kretschman = Riemann_contra[α,β,μ,ν] * Riemann_covar[α,β,μ,ν] 
 

        Kretschman_formula = Kretschmann_scalar(r,θ,a)
        @test isapprox(Kretschman,Kretschman_formula)

    end

end 






