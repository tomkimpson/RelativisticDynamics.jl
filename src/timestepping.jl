
#specify jacobian 

#https://www.youtube.com/watch?v=PxSALflWcE0
#ForwardDiff.jacobian - automatic differentiaon 
#solve with StaticArrays - whoch can work well for small system_parameters (< 20 ODEs)

function timestepping(X::PrognosticVariables, M::Model)

@unpack a, Tint = M.parameters
@unpack L,Q,m0,ϵ = M.constants

# Integration time 
tspan = (0.0,Tint) 


println("The model used is:")
println(M.parameters.model)

# Define the ODE to be solved
if M.parameters.model == :SphericalPhoton      
    params = [L,a]    
    u = vcat(X.xvector,X.pvector)             
    ode_prob = DifferentialEquations.ODEProblem(spherical_photon_hamiltonian!,u,tspan,params,progress = true)
    ode_solution = DifferentialEquations.solve(ode_prob,abstol=1e-8,reltol=1e-4)

elseif M.parameters.model == :MPD
    params = [a,m0,ϵ]
    u = vcat(X.xvector,X.pvector,X.svector)
    ode_prob = DifferentialEquations.ODEProblem(MPD!,u,tspan,params,progress = true)
    ode_solution = DifferentialEquations.solve(ode_prob,abstol=1e-8,reltol=1e-4)
end 





# Solve it 
#algorithm = DifferentialEquations.RK4() # probably define this elsewhere 
#ode_solution = DifferentialEquations.solve(ode_prob,algorithm,saveat=1)
#ode_solution = DifferentialEquations.solve(ode_prob,abstol=1e-8,reltol=1e-4)
#ode_solution=1.0
return ode_solution

    
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




function MPD!(du,u,p,τ)



    #Extract - can we do this better?
    t,r,θ,ϕ,pᵗ,pʳ,pᶿ,pᵠ,sᵗ,sʳ,sᶿ,sᵠ = u
    a,m0,ϵ = p
    xvector = [t,r,θ,ϕ]
    pvector = [pᵗ,pʳ,pᶿ,pᵠ]
    svector = [sᵗ,sʳ,sᶿ,sᵠ]


    # Define useful functions 
    g = covariant_metric(xvector,a)
    Σ = sigma(r,θ,a)
    Δ = delta(r,a)
    Γ = christoffel(r,θ,a)
    Riemann = riemann(r,θ,a) # This is the mixed contra/covar term
    metric_trace =-sin(θ)^2*Σ^2
    permutation_tensor = ϵ/sqrt(abs(metric_trace)) #This is also defined in the spintensor function 
    
    # Turn the contravariant permutation_tensor to eps^{ab}_{cd}
    @tensor begin
        eps[ρ,σ,μ,ν] := g[μ,λ]*g[ν,γ]*permutation_tensor[ρ,σ,λ,γ] 
    end

    @tensor begin
        Riemann_covar[μ,ν,ρ,σ] := g[μ,λ]*Riemann[λ,ν,ρ,σ] #This is the fully covariant form
    end

 
    Stensor = spintensor(xvector,pvector,svector,a,m0,ϵ)


    # 4 velocity 
    dx = calculate_four_velocity(pvector,Stensor,Riemann_covar,g,m0)


    # 4 momentum 
    dp = calculate_four_momentum(pvector,dx,svector,Γ,Riemann,eps,m0)
#

    #4- spin
    ds = calculate_four_spin(pvector,dx,svector,Γ,Riemann_covar,eps,m0)


    println("updating dx ds dp ")
    du[1:4] = dx
    du[5:8] = dp
    du[9:12]= ds

    # #Position 
    # du[1] = 0.0
    # du[2] = 0.0
    # du[3] = 0.0
    # du[4] = 0.0

    # #Momentum 
    # du[5] = 0.0
    # du[6] = 0.0
    # du[7] = 0.0
    # du[8] = 0.0


    # #Spin 
    # du[9]  = 0.0
    # du[10] = 0.0
    # du[11] = 0.0
    # du[12] = 0.0


    nothing #function returns nothing


end 




function calculate_four_velocity(pvector,Stensor,Riemann,g,m0)

    @tensor begin
        division_scalar = Riemann[μ,ν,ρ,σ]*Stensor[μ,ν]*Stensor[ρ,σ]
    end 
    
    @tensor begin
        dx[α] := -(pvector[α]+0.50*Stensor[α,β]*Riemann[β,γ,μ,ν]*pvector[γ]*Stensor[μ,ν]/(m0^2 + division_scalar))/m0^2
    end
    
    @tensor begin 
        Vsq = g[μ,ν]*dx[μ]*dx[ν]
    end 
    
    
    PV = -sqrt(-1.0/Vsq)
    dx = dx * PV
    
    
    #Check the value if needed
    @tensor begin 
        check_val = g[μ,ν]*dx[μ]*dx[ν]
    end 
    println("Check value")
    println(check_val)


    return dx
end 


function calculate_four_momentum(pvector,uvector,svector,Γ,Riemann,eps,m0)
     @tensor begin
        dp[α] := -Γ[α,μ,ν]*pvector[μ]*uvector[ν] + (Riemann[α,β,ρ,σ]*eps[ρ,σ,μ,ν]*svector[μ]*pvector[ν]*uvector[β])/(2.0*m0)
    end

    return dp
end 


function calculate_four_spin(pvector,uvector,svector,Γ,Riemann_covar,eps,m0)
    @tensor begin
       ds[α] := -Γ[α,μ,ν]*svector[μ]*uvector[ν] + (Riemann_covar[γ,β,ρ,σ]*eps[ρ,σ,μ,ν]*svector[μ]*pvector[ν]*svector[γ]*uvector[β])*pvector[α]/(2.0*m0^3)
   end

   return ds
end 