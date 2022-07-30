






#specify jacobian 

#https://www.youtube.com/watch?v=PxSALflWcE0
#ForwardDiff.jacobian - automatic differentiaon 
#solve with StaticArrays - whoch can work well for small system_parameters (< 20 ODEs)



function timestepping(X::PrognosticVariables, M::Model)


println("you are inside timestepping")


#u = [X.xvector[1],X.xvector[2],X.xvector[3],X.xvector[4],X.pvector[1],X.pvector[2],X.pvector[3],X.pvector[4],X.svector[1],X.svector[2],X.svector[3],X.svector[4]] # DEV: DONT DEFINE THIS AT EACH TIME STEP 

u = [1.0,2.0]
tspan = (0.0,1e5) # probably define this elsewhere
println(tspan)

ode_prob = DifferentialEquations.ODEProblem(MPD!,u,tspan)

algorithm = DifferentialEquations.RK4()
ode_solution = DifferentialEquations.solve(ode_prob,algorithm)

return ode_solution



# u0 = [1.0,0.0,0.0]
# tspan = (0.0,1e5)



# ode_prob = DifferentialEquations.ODEProblem(rober!,u0,tspan,params)
#ode_solve = DifferentialEquations.solve(ode_prob)
#println("Solution is:")
#println(ode_solve)
    
end 


# function rober!(du,u,p,t) #. This is an inplace function derivatives of the states, states, parameters, time 

# y1,y2,y3 = u
# k1,k2,k3 = p 
# du[1] = -k1*y1 + k3*y2*y3 
# du[2] = k1*y1 -k2*y2-k3*y2*y3
# du[3] = k2 * y2^2 

# nothing #function returns nothing

# end 


function MPD!(du,u,p,λ)

   # t,r,θ,ϕ,pᵗ,pʳ,pᶿ,pᵠ,sᵗ,sʳ,sᶿ,sᵠ = u 

    
    du[1] = 10.0
    du[2] = 20.0


    #Position 
    # du[1] = 0.0
    # du[2] = 0.0
    # du[3] = r^-1.50
    # du[4] = 0.0

    # #Momentum 
    # du[5] = 0.0
    # du[6] = 0.0
    # du[7] = 0.0
    # du[8] = 0.0

    # #Spin 
    # du[9] = 0.0
    # du[10] = 0.0
    # du[11] = 0.0
    # du[12] = 0.0
    

    nothing #function returns nothing


end 