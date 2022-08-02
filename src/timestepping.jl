






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
params = [Φ,Q,a,u0,u1]

u = vcat(X.xvector,X.pvector, X.svector)


# Define the ODE to be solved
ode_prob = DifferentialEquations.ODEProblem(MPD!,u,tspan,params)

# Solve it 
algorithm = DifferentialEquations.RK4() # probably define this elsewhere 
ode_solution = DifferentialEquations.solve(ode_prob,algorithm,saveat=1)

return ode_solution

    
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