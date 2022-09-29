
#specify jacobian 

#https://www.youtube.com/watch?v=PxSALflWcE0
#ForwardDiff.jacobian - automatic differentiaon 
#solve with StaticArrays - whoch can work well for small system_parameters (< 20 ODEs)

function timestepping(X::PrognosticVariables, M::Model)

@unpack a = M.parameters
@unpack L,Q,m0,ϵ,Tint = M.constants




println("The model used is:")
println(M.parameters.model)

println("Initial conditions are")
println(X.xvector)
println(X.pvector)
println(X.svector)

# Integration time 
tspan = (0.0,M.constants.Tint) 


#Bring all vectors together
u = vcat(X.xvector,X.pvector,X.svector)
params = [a,m0]
f = geodesic!

ode_prob = DifferentialEquations.ODEProblem(f,u,tspan,params)
ode_solution = DifferentialEquations.solve(ode_prob,DifferentialEquations.RK4()) # abstol=1e-9,reltol=1e-9 ,saveat = 1.0

return ode_solution
# # Define the ODE to be solved
# if M.parameters.model == :SphericalPhoton      
#     params = [L,a]    
#     u = vcat(X.xvector,X.pvector)
#     f = spherical_photon_hamiltonian!     
# elseif M.parameters.model == :MPD
#     params = [a,m0,ϵ]
#     u = vcat(X.xvector,X.pvector,X.svector)
#     #f = MPD!
#     f = geodesic!
# end 

# println("INTEGRATING!")
# println("ODE PROB")
# ode_prob = DifferentialEquations.ODEProblem(f,u,tspan,params,progress = true) #progress_steps = 1
# println("ODE SOLN")
# ode_solution = DifferentialEquations.solve(ode_prob,DifferentialEquations.RK4()) # abstol=1e-9,reltol=1e-9 ,saveat = 1.0
# return ode_solution

    
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

    println(r,"  ", uvector[2])
    

    nothing #function returns nothing


end 



















# function derivatives(y,parameters)


#     a,m0 = parameters

    


#     dy = zeros(Float64,8) 

#     return dy 

# end 















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


# function geodesic!(du,u,p,λ)

#     t,r,θ,ϕ,pᵗ,pʳ,pᶿ,pᵠ = u
#     L,a = p


#     # Define useful functions 
#     Σ = sigma(r,θ,a)
#     Δ = delta(r,a)
#     κ = pᶿ^2 +L^2/sin(θ)^2 + a^2*sin(θ)^2

#     #Position 
#     du[1] = 1.0+(2.0*r*(r^2 +a^2)-2*a*r*L)/(Σ*Δ)
#     du[2] = Δ*pʳ/Σ
#     du[3] = pᶿ/Σ
#     du[4] = (2.0*r*a + (Σ - 2.0*r)*L/sin(θ)^2) / (Σ*Δ)

#     #Momentum 
#     du[5] = 0.0
#     du[6] = (-κ*(r-1) + 2*r*(r^2+a^2) - 2*a*L)/(Σ*Δ) - 2*pʳ^2*(r-1)/Σ
#     du[7] = sin(θ)*cos(θ)*(L^2/sin(θ)^4 -a^2)/Σ
#     du[8] = 0.0

#     nothing #function returns nothing


# end 



















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









# ScRATCH spACE

# #t
# dp1_a = -Γ[1,1,2]*pvector[1]*velocity[2] - Γ[1,2,1]*pvector[2]*velocity[1]
# dp1_b = -Γ[1,1,3]*pvector[1]*velocity[3] - Γ[1,3,1]*pvector[3]*velocity[1]
# dp1_c = -Γ[1,2,4]*pvector[2]*velocity[4] - Γ[1,4,2]*pvector[4]*velocity[2]
# dp1_d = -Γ[1,3,4]*pvector[3]*velocity[4] - Γ[1,4,3]*pvector[4]*velocity[3]
# dp1 = dp1_a+dp1_b+dp1_c+dp1_d 

# #r
# dp2_a = -Γ[2,1,1]*pvector[1]*velocity[1] 
# dp2_b = -Γ[2,1,4]*pvector[1]*velocity[4] -Γ[2,4,1]*pvector[4]*velocity[1]
# dp2_c = -Γ[2,2,2]*pvector[2]*velocity[2] 
# dp2_d = -Γ[2,2,3]*pvector[2]*velocity[3] -Γ[2,3,2]*pvector[3]*velocity[2]
# dp2_e = -Γ[2,3,3]*pvector[3]*velocity[3] 
# dp2_f = -Γ[2,4,4]*pvector[4]*velocity[4]
# dp2 = dp2_a+dp2_b+dp2_c+dp2_d + dp2_e+dp2_f 

# #θ
# dp3_a = -Γ[3,1,1]*pvector[1]*velocity[1] 
# dp3_b = -Γ[3,1,4]*pvector[1]*velocity[4] -Γ[3,4,1]*pvector[4]*velocity[1]
# dp3_c = -Γ[3,2,2]*pvector[2]*velocity[2] 
# dp3_d = -Γ[3,2,3]*pvector[2]*velocity[3] -Γ[3,3,2]*pvector[3]*velocity[2]
# dp3_e = -Γ[3,3,3]*pvector[3]*velocity[3] 
# dp3_f = -Γ[3,4,4]*pvector[4]*velocity[4]
# dp3 = dp3_a+dp3_b+dp3_c+dp3_d + dp3_e+dp3_f 

# #ϕ
# dp4_a = -Γ[4,1,2]*pvector[1]*velocity[2] - Γ[4,2,1]*pvector[2]*velocity[1]
# dp4_b = -Γ[4,1,3]*pvector[1]*velocity[3] - Γ[4,3,1]*pvector[3]*velocity[1]
# dp4_c = -Γ[4,3,4]*pvector[3]*velocity[4] - Γ[4,4,3]*pvector[4]*velocity[3]
# dp4_d = -Γ[4,2,4]*pvector[2]*velocity[4] - Γ[4,4,2]*pvector[4]*velocity[2]
# dp4 = dp4_a+dp4_b+dp4_c+dp4_d 
