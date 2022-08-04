






#specify jacobian 

#https://www.youtube.com/watch?v=PxSALflWcE0
#ForwardDiff.jacobian - automatic differentiaon 
#solve with StaticArrays - whoch can work well for small system_parameters (< 20 ODEs)



function timestepping(X::PrognosticVariables, M::Model)



println("dev: you are inside timestepping")
@unpack a, Tint = M.parameters
@unpack Φ,Q,u0,u1 = M.constants


# Integration time 
tspan = (0.0,Tint) 

#Params
# params = [Φ,Q,a,u0,u1]
# u = vcat(X.xvector,X.pvector, X.svector)



u = vcat(X.xvector,X.pvector)
params = [Φ,a]

println("dev: this is uvec", u)
println("dev: this is params", params)


# Define the ODE to be solved
ode_prob = DifferentialEquations.ODEProblem(spherical_photon_hamiltonian!,u,tspan,params)

# Solve it 
#algorithm = DifferentialEquations.RK4() # probably define this elsewhere 
#ode_solution = DifferentialEquations.solve(ode_prob,algorithm,saveat=1)
ode_solution = DifferentialEquations.solve(ode_prob,abstol=1e-8,reltol=1e-4)

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
















function MPD!(du,u,p,λ)

    t,r,χ,ϕ,pᵗ,pʳ,pᶿ,pᵠ,sᵗ,sʳ,sᶿ,sᵠ = u     
    Φ,Q,a,u0,u1 = p

    # Define θ from χ
    θ = acos(sqrt(u0) * sin(χ))
    

    # Define useful functions 
    Σ = sigma(r,θ,a)
    Δ = delta(r,a)

    #Position 
    du[1] = 0.0
    du[2] = 0.0
    du[3] = a*sqrt(u0 * sin(χ)^2 - u1) / (r^2 + a^2 * u0 * sin(χ)^2)
    du[4] = (2.0*r*a + (Σ - 2.0*r)*Φ/sin(θ)^2) / (Σ*Δ)

    #Momentum 
    du[5] = 0.0
    du[6] = 0.0
    du[7] = 0.0
    du[8] = 0.0

    #Spin 
    du[9] = 0.0
    du[10] = 0.0
    du[11] = 0.0
    du[12] = 0.0
    

    nothing #function returns nothing


end 