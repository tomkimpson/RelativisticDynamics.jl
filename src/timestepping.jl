
using DifferentialEquations
using ComponentArrays
using SciMLSensitivity
using Parameters: @unpack


"""
The timesteppig Integration
"""
function timestepping(X::PrognosticVariables, M::Model)

@unpack a,e = M.parameters
@unpack m0, Tint = M.constants
# @unpack L,Q,m0,Tint = M.constants


#tstart:: Real = 0.0
#tend::
# # Integration time 
tspan = (0.0,M.constants.Tint) 

println("tint type:", typeof(Tint))
# #Bring all vectors together
u = vcat(X.xvector,X.pvector,X.svector)

params = [a,m0]

println("TYPES - u1 ",typeof(u[1]))
println("TYPES - u ",typeof(u))
println("TYPES - t ",typeof(tspan))
println("TYPES - p ",typeof(params))



f = MPD! #The ODE 
println("HERE")
ode_prob = ODEProblem(MPD!,u,tspan,params)
ode_solution = solve(ode_prob,DifferentialEquations.RK4()) # abstol=1e-9,reltol=1e-9 ,saveat = 1.0

return ode_solution



#p = X.pvector*1e7; u0 = [1.0;1.0]
#prob = ODEProblem(lotka_volterra!,u0,(0.0,10.0),p)
#sol = solve(prob,Tsit5(),reltol=1e-6,abstol=1e-6)

end 



function lotka_volterra!(du,u,p,t)
  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] = dy = -p[3]*u[2] + p[4]*u[1]*u[2]
end





function MPD!(du,u,p,τ)

    #println("isnide MPD")
    #Extract the coordinates/constants 
    t,r,θ,ϕ,pᵗ,pʳ,pᶿ,pᵠ,sᵗ,sʳ,sᶿ,sᵠ = u # coordinate variables
    a,m0 = p                            # constants

    #Define vectors from the coordinate variables. 
    xvector = [t,r,θ,ϕ]
    pvector = [pᵗ,pʳ,pᶿ,pᵠ]
    svector = [sᵗ,sʳ,sᶿ,sᵠ]


    # # #Define some useful quantities for this timestep 
    g = covariant_metric(xvector,a)        # the metric 
    Γ = christoffel(xvector,a)            # the Christoffel symbols
    Riemann = riemann(xvector,a)          #the mixed contra/covar Riemann term R^{a}_{bcd}
    @tullio Riemann_covar[μ,ν,ρ,σ] := g[μ,λ]*Riemann[λ,ν,ρ,σ] #This is the fully covariant form R_{abcd}
    



    levi = permutation_tensor(g)  #This is the fully contravariant Levi Civita tensor 
    @tullio levi_mixed[ρ,σ,μ,ν] := g[μ,x]*g[ν,y] * levi[ρ,σ,x,y]
    
  
    stensor = spintensor(levi,pvector,svector,m0) #the fully contravariant spin tensor s^{ab}


    # #Get the derivative objects
    # #4-velocity
    uvector = calculate_four_velocity(pvector,stensor,Riemann_covar,g,m0)

    # #4-momentum
    dp = calculate_four_momentum(pvector,uvector,svector,Γ,Riemann,levi_mixed,m0)

    # #4-spin
    ds = calculate_four_spin(pvector,uvector,svector,Γ,Riemann_covar,levi_mixed,m0)

    #println(uvector)
    #Package it up 
    du[1:4] = uvector
    du[5:8] = dp
    du[9:12]= ds
    # println("outputs")
    # println(uvector)
    # println(dp)
    # println(ds)
    # println("-----------------------------------")


    #du[1:4] = [0.0,pʳ*g[1,1],0.0,0.0]
    #du[5:8] = [0.0,0.0,0.0,0.0]
    #du[9:12]= [0.0,0.0,0.0,0.0]

    nothing #function returns nothing


end 


function calculate_four_velocity(pvector,Stensor,Riemann,g,m0)
    #println("enter 4 velocity")
    #println(g)

    @tullio scalar_divisor := Riemann[μ,ν,ρ,σ]*Stensor[μ,ν]*Stensor[ρ,σ] / 4.0
    

    @tullio correction[α] := 0.50 *(Stensor[α,β]*Riemann[β,γ,μ,ν]*pvector[γ]*Stensor[μ, ν])/(m0^2 + scalar_divisor)
    

    
    @tullio dx[α] := -(pvector[α] + correction[α])/m0^2 
    #@tullio dx[α] := -(pvector[α] )/m0^2 

    

    @tullio  Vsq := g[μ,ν]*dx[μ]*dx[ν] 

#    println("VSQ = ")
#    println(Vsq)
#    println("dx = ")
#    println(dx)
    PV = -sqrt(-1.0/Vsq)
    dx = dx * PV

    return dx
end 


function calculate_four_momentum(pvector,uvector,svector,Γ,Riemann,levi_mixed,m0)



    @tullio correction[α] := (Riemann[α,β,ρ,σ]*levi_mixed[ρ,σ,μ,ν]*svector[μ]*pvector[ν]*uvector[β])/(2*m0)



    @tullio dp[α] := -Γ[α,μ,ν]*pvector[μ]*uvector[ν] + correction[α]
   

   return dp
end 


function calculate_four_spin(pvector,uvector,svector,Γ,Riemann_covar,levi_mixed,m0)


    @tullio correction[α] := pvector[α]*(Riemann_covar[γ,β,ρ,σ]*levi_mixed[ρ,σ,μ,ν]*svector[μ]*pvector[ν]*svector[γ]*uvector[β])/(2*m0^3)
    

    @tullio ds[α] := -Γ[α,μ,ν]*svector[μ]*uvector[ν] + correction[α]
   

   return ds
end 



