






#specify jacobian 

#https://www.youtube.com/watch?v=PxSALflWcE0
#ForwardDiff.jacobian - automatic differentiaon 
#solve with StaticArrays - whoch can work well for small system_parameters (< 20 ODEs)



function timestepping(X::PrognosticVariables, M::Model)

@unpack a, Tint = M.parameters
@unpack Φ,Q,u0,u1,m0,ϵ = M.constants

# Integration time 
tspan = (0.0,Tint) 


println("The model used is:")
println(M.parameters.model)

# Define the ODE to be solved
if M.parameters.model == :SphericalPhoton      
    params = [Φ,a]    
    u = vcat(X.xvector,X.pvector)             
    ode_prob = DifferentialEquations.ODEProblem(spherical_photon_hamiltonian!,u,tspan,params,progress = true)
    ode_solution = DifferentialEquations.solve(ode_prob,abstol=1e-8,reltol=1e-4)

elseif M.parameters.model == :MPD
    params = [Φ,a,m0,ϵ]
    u = vcat(X.xvector,X.pvector,X.svector)
    #ode_prob = DifferentialEquations.ODEProblem(MPD!,u,tspan,params,progress = true)
    #ode_solution = DifferentialEquations.solve(ode_prob,abstol=1e-8,reltol=1e-4)

    #MPD!(du,u,params,τ)
    println("MPD TEST")
    MPD_test(u,params)

end 








# Solve it 
#algorithm = DifferentialEquations.RK4() # probably define this elsewhere 
#ode_solution = DifferentialEquations.solve(ode_prob,algorithm,saveat=1)
#ode_solution = DifferentialEquations.solve(ode_prob,abstol=1e-8,reltol=1e-4)
ode_solution=1.0
return ode_solution

    
end 



function spherical_photon_hamiltonian!(du,u,p,λ)

    t,r,θ,ϕ,pᵗ,pʳ,pᶿ,pᵠ = u
    Φ,a = p


    # Define useful functions 
    Σ = sigma(r,θ,a)
    Δ = delta(r,a)

    #Position 
    du[1] = 0.0
    du[2] = 0.0
    du[3] = pᶿ/Σ
    du[4] = (2.0*r*a + (Σ - 2.0*r)*Φ/sin(θ)^2) / (Σ*Δ)

    #Momentum 
    du[5] = 0.0
    du[6] = 0.0
    du[7] = sin(θ)*cos(θ)*(Φ^2/sin(θ)^4 -a^2)/Σ
    du[8] = 0.0

    nothing #function returns nothing


end 




function MPD!(du,u,p,τ)

    #Extract - can we do this better?
    t,r,θ,ϕ,pᵗ,pʳ,pᶿ,pᵠ,sᵗ,sʳ,sᶿ,sᵠ = u
    Φ,a,m0,ϵ = p
    xvector = [t,r,θ,ϕ]
    pvector = [pᵗ,pʳ,pᶿ,pᵠ]
    svector = [sᵗ,sʳ,sᶿ,sᵠ]


    # Define useful functions 
    g = covariant_metric(xvector,a)
    Σ = sigma(r,θ,a)
    Δ = delta(r,a)
    Γ = christoffel(r,θ,a)
    Riemann = riemann(r,θ,a)

 
    Stensor = spintensor(xvector,pvector,svector,a,m0,ϵ)


    # 4 velocity 
    @tensor begin
    division_scalar = Riemann[μ,ν,ρ,σ]*Stensor[μ,ν]*Stensor[ρ,σ]
    end 


    



    @tensor begin
        dx[α] := -(pvector[α]+0.50*Stensor[α,β]*Riemann[β,γ,μ,ν]*pvector[γ]*Stensor[μ,ν]/(m0^2 + division_scalar))/m0^2
    end



    @tensor begin 
        Vsq = g[μ,ν]*dx[μ]*dx[ν]
    end 

    PV = -sqrt(1.0/Vsq)

    println("Vsq:")
    println(Vsq)
    println("PV:")
    println(PV)


    #dx = dx * PV
    display(dx)
    println("CHECK")
    @tensor begin
        check_val = g[μ,ν]*dx[μ]*dx[ν]
    end 
    println(check_val)    


    return



    #Position 
    du[1] = 0.0
    du[2] = 0.0
    du[3] = pᶿ/Σ
    du[4] = (2.0*r*a + (Σ - 2.0*r)*Φ/sin(θ)^2) / (Σ*Δ)

    #Momentum 
    du[5] = 0.0
    du[6] = 0.0
    du[7] = sin(θ)*cos(θ)*(Φ^2/sin(θ)^4 -a^2)/Σ
    du[8] = 0.0

    nothing #function returns nothing


end 


function MPD_test(u,p)

    #Extract - can we do this better?
    t,r,θ,ϕ,pᵗ,pʳ,pᶿ,pᵠ,sᵗ,sʳ,sᶿ,sᵠ = u
    Φ,a,m0,ϵ = p
    xvector = [t,r,θ,ϕ]
    pvector = [pᵗ,pʳ,pᶿ,pᵠ]
    svector = [sᵗ,sʳ,sᶿ,sᵠ]


    # Define useful functions 
    g = covariant_metric(xvector,a)
    Σ = sigma(r,θ,a)
    Δ = delta(r,a)
    Γ = christoffel(r,θ,a)
    Riemann = riemann(r,θ,a)

 
    Stensor = spintensor(xvector,pvector,svector,a,m0,ϵ)


    # 4 velocity 
    @tensor begin
    division_scalar = Riemann[μ,ν,ρ,σ]*Stensor[μ,ν]*Stensor[ρ,σ]
    end 


    



    @tensor begin
        dx[α] := -(pvector[α]+0.50*Stensor[α,β]*Riemann[β,γ,μ,ν]*pvector[γ]*Stensor[μ,ν]/(m0^2 + division_scalar))/m0^2
    end



    @tensor begin 
        Vsq = g[μ,ν]*dx[μ]*dx[ν]
    end 

    PV = -sqrt(1.0/Vsq)

    println("pvector")
    println(pvector)

    println("Vsq:")
    println(Vsq)
    println("PV:")
    println(PV)


    #dx = dx * PV
    display(dx)
    println("CHECK")
    @tensor begin
        check_val = g[μ,ν]*dx[μ]*dx[ν]
    end 
    println(check_val)    


    nothing #function returns nothing


end 


