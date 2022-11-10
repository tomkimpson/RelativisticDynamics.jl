
using DifferentialEquations
using ComponentArrays
using SciMLSensitivity
using Parameters: @unpack


"""
    solution = timestepping(X,M)
The timestepping integration once all variables have been initialised
"""
function timestepping(X::PrognosticVariables, M::Model)

@unpack a,e = M.parameters
@unpack m0, Tint, integrator_method = M.constants


tspan = (0.0,M.constants.Tint) 
u = vcat(X.xvector,X.pvector,X.svector)
params = [a,m0]

f = MPD! #The ODE 
ode_prob = ODEProblem(MPD!,u,tspan,params)
ode_solution = solve(ode_prob,DifferentialEquations.RK4())



#h = 1e-15

#τ = 0.0 
#println(M.constants.Tint)

# y = u
# IO_array = zeros(Float64,trunc(Int,2.0*M.constants.Tint/0.5),12)
# i = 1

# IO_array[i,:] = y
# i += 1


# while τ < M.constants.Tint

#     ynew,h,τ,saveit = RK_step(y,h,τ,params)

#     if saveit 
#         y = ynew
#         IO_array[i,:] = y
#         i += 1



#     end 
#     #println(ynew[2], " ", τ," ",M.constants.Tint," ", h)


# end 

# rows_to_keep = findall(IO_array[:,1] .!= 0)
# println(rows_to_keep)
# ode_solution = IO_array[rows_to_keep,:]
return ode_solution


end 




function RK_step(u_in,h,τ,p)


    B21 =1.0/5.0
    B31 = 3.0/40.0
    B32 = 9.0/40.0
    B41 = 3.0/10.0
    B42 = -9.0/10.0
    B43 = 6.0/5.0 
    B51 = -11.0/54.0
    B52 = 5.0/2.0
    B53 = -70.0/27.0
    B54 = 35.0/27.0
    B61 = 1631.0/55296.0
    B62 = 175.0/512.0
    B63 = 575.0/13824.0
    B64 = 44275.0/110592.0
    B65 = 253.0/4096.0
    c1 = 37.0/378.0
    c3 = 250.0/621.0
    c4 = 125.0/594.0
    c6=512.0/1771.0
    cbar1 = 2825.0/27648.0
    cbar3 = 18575.0/48384.0
    cbar4=13525.0/55296.0
    cbar5 = 277.0/14336.0
    cbar6 = 1.0/4.0

    y1 = u_in 
    du = similar(y1)
    MPD!(du,y1,p,τ)
    k1 = h * du
   


    y2 = y1 + B21*k1 
    MPD!(du,y2,p,τ)
    k2 = h * du


    y3 = y1 + B31 *k1 + B32*k2 
    MPD!(du,y3,p,τ)
    k3 = h * du


    y4 = y1 + B41 *k1 + B42*k2 + B43*k3 
    MPD!(du,y4,p,τ)
    k4 = h * du

    y5 = y1 + B51 *k1 + B52*k2 + B53*k3 + B54*k4 
    MPD!(du,y5,p,τ)
    k5 = h * du
    

    y6 = y1 + B61 *k1 + B62*k2 + B63*k3 + B64*k4 + B65*k5 
    MPD!(du,y6,p,τ)
    k6 = h * du

   
    ynew = y1 + c1*k1  + c3*k3 + c4*k4  +c6*k6 
    yerr = y1 + cbar1*k1 + cbar3*k3 + cbar4*k4 + cbar5*k5 + cbar6*k6


    Δ = abs.(ynew - yerr)
    yscal = abs.(y1) + abs.(k1) 

    ratio = Δ ./ yscal

    escal = 1e15

    errmax = escal * maximum(ratio)
    # println(escal)
    # println(maximum(ratio))
    # println(ratio)
    # println("---------------------")
 
    if errmax > 1.0 
        #println("errmax > 1")
        h = ShrinkStepsize(errmax,h)
        return y1, h,τ,false
    else
        
        τ += h
        h = GrowStepsize(errmax,h)
        #println("errmax < 1")
        #println(h)
        return ynew,h,τ,true 
    end 

end 



function GrowStepsize(errmax,h)


    S = 0.90
    Pgrow = -0.20
    Pshrink = -0.250
    errcon = (5.0/S)^(1.0/Pgrow)


    if errmax > errcon 
        h = S*h*errmax^Pgrow
    else
        h = h * 5.0
    end

    return h 
end 





function ShrinkStepsize(errmax,h)

    S = 0.90
    Pshrink = -0.250
    htemp = S*h*errmax^Pshrink
    return copysign(max(abs(htemp),0.10*abs(h)),h)

end 










function MPD!(dx,u,p,τ)

    # println("MPD FUNC IN")
    # println(u)
    # println("--------------------")

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
    #println(u)
    #println(dx)
    #Package it up 
    dx[1:4] = uvector
    #println(u)
    dx[5:8] = dp
    dx[9:12]= ds
    #println(u)

    # println("MPD FUNC OUT")
    # println(u)
    # println(dx)
    # println("--------------------")
 
    nothing #function returns nothing


end 


function calculate_four_velocity(pvector,Stensor,Riemann,g,m0)

    @tullio scalar_divisor := Riemann[μ,ν,ρ,σ]*Stensor[μ,ν]*Stensor[ρ,σ] / 4.0
    
    @tullio correction[α] := 0.50 *(Stensor[α,β]*Riemann[β,γ,μ,ν]*pvector[γ]*Stensor[μ, ν])/(m0^2 + scalar_divisor)
        
    @tullio dx[α] := -(pvector[α] + correction[α])/m0^2 

    @tullio  Vsq := g[μ,ν]*dx[μ]*dx[ν] 

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



