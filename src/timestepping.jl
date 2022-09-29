
#specify jacobian 

#https://www.youtube.com/watch?v=PxSALflWcE0
#ForwardDiff.jacobian - automatic differentiaon 
#solve with StaticArrays - whoch can work well for small system_parameters (< 20 ODEs)

function timestepping(X::PrognosticVariables, M::Model)

@unpack a = M.parameters
@unpack L,Q,m0,ϵ,Tint = M.constants


# Integration time 
tspan = (0.0,M.constants.Tint) 


#Bring all vectors together
u = vcat(X.xvector,X.pvector,X.svector)
params = [a,m0]
f = geodesic!

ode_prob = DifferentialEquations.ODEProblem(f,u,tspan,params)
ode_solution = DifferentialEquations.solve(ode_prob,DifferentialEquations.RK4()) # abstol=1e-9,reltol=1e-9 ,saveat = 1.0

return ode_solution
    
end 



function geodesic!(du,u,p,τ)

    #Extract - can we do this better?
    t,r,θ,ϕ,pᵗ,pʳ,pᶿ,pᵠ,sᵗ,sʳ,sᶿ,sᵠ = u #variables
    a,m0 = p                          #parameters
    
    #Define 4 momentum and 4 velocity vectors
    pvector = [pᵗ,pʳ,pᶿ,pᵠ]
    uvector = pvector/m0

    #What is the EOM for the 4-momentum derivative?
    Γ = christoffel(r,θ,a)
    @tensor begin
        dp[α] := -Γ[α,μ,ν]*pvector[μ]*uvector[ν] 
    end


    du[1:4] = uvector
    du[5:8] = dp
    du[9:12]= [0.0,0.0,0.0,0.0]

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


    # println("dx = ",dx)
    # println("pv = ",pvector)
    # println("dp = ",dp)
    # println("ds = ",ds)
    # println("----------------------------")

    #println("updating dx ds dp ")
    du[1:4] = dx
    du[5:8] = dp
    du[9:12]= ds


    #println(du)

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
        dx[α] := -(pvector[α])/m0^2 #+0.50*Stensor[α,β]*Riemann[β,γ,μ,ν]*pvector[γ]*Stensor[μ,ν]/(m0^2 + division_scalar))/m0^2
    end
    
    @tensor begin 
        Vsq = g[μ,ν]*dx[μ]*dx[ν]
    end 
    
    
   # PV = -sqrt(-1.0/Vsq)
   # dx = dx * PV
    
    
    #Check the value if needed
    # @tensor begin 
    #     check_val = g[μ,ν]*dx[μ]*dx[ν]
    # end 
    # println("Check value")
    # println(check_val)


    return dx
end 


function calculate_four_momentum(pvector,uvector,svector,Γ,Riemann,eps,m0)
     @tensor begin
        dp[α] := -Γ[α,μ,ν]*pvector[μ]*uvector[ν] #+ (Riemann[α,β,ρ,σ]*eps[ρ,σ,μ,ν]*svector[μ]*pvector[ν]*uvector[β])/(2.0*m0)
    end

    return dp
end 


function calculate_four_spin(pvector,uvector,svector,Γ,Riemann_covar,eps,m0)
    @tensor begin
       ds[α] := -Γ[α,μ,ν]*svector[μ]*uvector[ν]# + (Riemann_covar[γ,β,ρ,σ]*eps[ρ,σ,μ,ν]*svector[μ]*pvector[ν]*svector[γ]*uvector[β])*pvector[α]/(2.0*m0^3)
   end

   return ds
end 









