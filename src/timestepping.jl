

"""
The timesteppig Integration
"""
function timestepping(X::PrognosticVariables, M::Model)

@unpack a = M.parameters
@unpack L,Q,m0,ϵ,Tint = M.constants


# Integration time 
tspan = (0.0,M.constants.Tint) 


#Bring all vectors together
u = vcat(X.xvector,X.pvector,X.svector)
params = [a,m0,ϵ]
f = MPD! #The ODE 

ode_prob = DifferentialEquations.ODEProblem(f,u,tspan,params)
ode_solution = DifferentialEquations.solve(ode_prob,DifferentialEquations.RK4()) # abstol=1e-9,reltol=1e-9 ,saveat = 1.0

return ode_solution
    
end 



function MPD!(du,u,p,τ)

    #Extract the coordinates/constants 
    t,r,θ,ϕ,pᵗ,pʳ,pᶿ,pᵠ,sᵗ,sʳ,sᶿ,sᵠ = u # coordinate variables
    a,m0,ϵ = p                          # constants 
    
    #Define vectors from the coordinate variables. 
    xvector = [t,r,θ,ϕ]
    pvector = [pᵗ,pʳ,pᶿ,pᵠ]
    svector = [sᵗ,sʳ,sᶿ,sᵠ]


    #Define some useful quantities for this timestep 
    g = covariant_metric(xvector,a)   #the metric 
    Γ = christoffel(r,θ,a)            #the Christoffel symbols
    Riemann = riemann(r,θ,a)          #the mixed contra/covar Riemann term R^{a}_{bcd}
    @tensor begin
        Riemann_covar[μ,ν,ρ,σ] := g[μ,λ]*Riemann[λ,ν,ρ,σ] #This is the fully covariant form R_{abcd}
    end

    levi = permutation_tensor(g,ϵ)  #This is the fully contravariant Levi Civita tensor 
    @tensor begin 
        levi_mixed[ρ,σ,μ,ν] := g[μ,x]*g[ν,y] * levi[ρ,σ,x,y]
    end 
  
    stensor = spintensor(levi,pvector,svector,m0) #the fully contravariant spin tensor s^{ab}


    #Get the derivative objects

    #4-velocity
    uvector = calculate_four_velocity(pvector,stensor,Riemann_covar,g,m0)

    #4-momentum
    dp = calculate_four_momentum(pvector,uvector,svector,Γ,Riemann,levi_mixed,m0)

    #4-spin
    ds = calculate_four_spin(pvector,uvector,svector,Γ,Riemann_covar,levi_mixed,m0)

    #Package it up 
    du[1:4] = uvector
    du[5:8] = dp
    du[9:12]= ds

    nothing #function returns nothing


end 




function calculate_four_velocity(pvector,Stensor,Riemann,g,m0)

    @tensor begin 
        scalar_divisor = Riemann[μ,ν,ρ,σ]*Stensor[μ,ν]*Stensor[ρ,σ] / 4.0
    end 

    @tensor begin
        correction[α] := 0.50 *(Stensor[α,β]*Riemann[β,γ,μ,ν]*pvector[γ]*Stensor[μ, ν])/(m0^2 + scalar_divisor)
    end

    
    @tensor begin
        dx[α] := -(pvector[α] + correction[α])/m0^2 
    end

    @tensor begin 
        Vsq = g[μ,ν]*dx[μ]*dx[ν]
    end 

    PV = -sqrt(-1.0/Vsq)
    dx = dx * PV
    
    return dx
end 


function calculate_four_momentum(pvector,uvector,svector,Γ,Riemann,levi_mixed,m0)



    @tensor begin
        correction[α] := (Riemann[α,β,ρ,σ]*levi_mixed[ρ,σ,μ,ν]*svector[μ]*pvector[ν]*uvector[β])/(2*m0)

    end


    @tensor begin
       dp[α] := -Γ[α,μ,ν]*pvector[μ]*uvector[ν] + correction[α]
   end

   return dp
end 


function calculate_four_spin(pvector,uvector,svector,Γ,Riemann_covar,levi_mixed,m0)


    @tensor begin
        correction[α] := pvector[α]*(Riemann_covar[γ,β,ρ,σ]*levi_mixed[ρ,σ,μ,ν]*svector[μ]*pvector[ν]*svector[γ]*uvector[β])/(2*m0^3)
    end



    @tensor begin
       ds[α] := -Γ[α,μ,ν]*svector[μ]*uvector[ν] + correction[α]
   end

   return ds
end 



