






#specify jacobian 

#https://www.youtube.com/watch?v=PxSALflWcE0
#ForwardDiff.jacobian - automatic differentiaon 
#solve with StaticArrays - whoch can work well for small system_parameters (< 20 ODEs)



function timestepping(X::PrognosticVariables, M::Model)



println("dev: you are inside timestepping")
@unpack a = M.parameters

#u = vcat(X.xvector, X.pvector,X.svector)


####-------- Spherical photon orbits 
# Pick a radius and an initial θ, ϕ
rmin = 2.0 * (1.0 + cos(2.0/3.0 * acos(-a)))
rmax = 2.0 * (1.0 + cos(2.0/3.0 * acos(a)))
r = (rmin+rmax)/2.0
θ = π/6.0
ϕ = 0.0

# Define constants of motion 
Σ = sigma(r,θ,a)
Φ = - (r^3 - 3.0*r^2 + a^2*r + a^2) / (a*(r - 1.0))
Q = -(r^3 * (r^3 - 6.0*r^2 + 9.0 * r - 4.0*a^2))/(a^2 * (r - 1.0)^2)


# Define u0 function 
u0 = ((a^2 - Q - Φ^2) + sqrt((a^2 - Q - Φ^2)^2 + 4.0*a^2 * Q))/(2.0*a^2)
u1 = ((a^2 - Q - Φ^2) - sqrt((a^2 - Q - Φ^2)^2 + 4.0*a^2 * Q))/(2.0*a^2)





tspan = (0.0,100) # probably define this elsewhere




params = [Φ,Q,a,u0,u1]
χ = acos(cos(θ)/sqrt(u0))
temp_vector = [0.0,r,χ,0.0]
u = vcat(temp_vector,X.pvector, X.svector)




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