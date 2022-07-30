






#specify jacobian 

#https://www.youtube.com/watch?v=PxSALflWcE0
#ForwardDiff.jacobian - automatic differentiaon 
#solve with StaticArrays - whoch can work well for small system_parameters (< 20 ODEs)



function timestepping(X::PrognosticVariables, M::Model)


println("dev: you are inside timestepping")
u = vcat(X.xvector, X.pvector,X.svector)
tspan = (0.0,100) # probably define this elsewhere




# Define the ODE to be solved
ode_prob = DifferentialEquations.ODEProblem(MPD!,u,tspan)

# Solve it 
algorithm = DifferentialEquations.RK4() # probably define this elsewhere 
ode_solution = DifferentialEquations.solve(ode_prob,algorithm,saveat=1)

return ode_solution

    
end 


function MPD!(du,u,p,λ)

    t,r,θ,ϕ,pᵗ,pʳ,pᶿ,pᵠ,sᵗ,sʳ,sᶿ,sᵠ = u 



    #Position 
    du[1] = 0.0
    du[2] = 0.0
    du[3] = 0.0
    du[4] = 10.0

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