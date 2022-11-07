
"""
The main output function from this package.
"""
function orbit(::Type{NF}=Float64;              # number format, use Float64 as default
    kwargs...                   # all additional non-default parameters
    ) where {NF<:AbstractFloat}




    println("NF IS =", NF)



    # Setup all system parameters, universal constants etc.
    P = SystemParameters(NF=NF;kwargs...) # Parameters
    bounds_checks(P)

    println(P.a, typeof(P.a))

    C = Constants(P)                      # Constants
    M = Model(P,C)                        # Pack all of the above into a single *Model struct 

    println("Called orbit() with e = ", P.e)


    #Initial conditions 
    initialization = initial_conditions(M)

    #Evolve in time
    solution = timestepping(initialization, M)
    
    return solution, M

end



# using Zygote, SciMLSensitivity
# using DifferentialEquations

# function lotka_volterra!(du,u,p,t)
#   du[1] = p[1]*u[1] - p[2]*u[1]*u[2]
#   du[2] = -p[3]*u[2] + p[4]*u[1]*u[2]
# end

# function MPD!(du,u,p,τ)

#     #Extract the coordinates/constants 
#     t,r,θ,ϕ,pᵗ,pʳ,pᶿ,pᵠ,sᵗ,sʳ,sᶿ,sᵠ = u 
#     a,m0 = p

#     du[1:4] = [0.0,pʳ,0.0,0.0]
#     du[5:8] = [0.0,0.0,0.0,0.0]
#     du[9:12]= [0.0,0.0,0.0,0.0]

# end 





# function orbit(u0,p,x)



#     P = SystemParameters(NF=Float64) # Parameters

#     #p = [1.5,1.0,3.0,1.0]; u0 = [1.0;1.0]
    
    
    
#     prob = ODEProblem(lotka_volterra!,u0,(0.0,10.0),p)
#     sum(solve(prob,Tsit5(),reltol=1e-6,abstol=1e-6,saveat=0.1))


# end 



"""
Look before you leap
"""
function bounds_checks(P::SystemParameters)
    @boundscheck 0.0<=P.a<1.0          || throw(error("Spin parameter a is out of bounds"))
    @boundscheck P.mBH/P.mPSR >= 1e3   || throw(error("Mass ratio is too small"))
    @boundscheck 1 <= P.rPSR <= 100    || throw(error("Pulsar radius is unphysical"))
    @boundscheck 0.0<  P.e<1.0         || throw(error("Eccentricity is outside range"))
    @boundscheck P.orbit_dir in [-1,1] || throw(error("Orbit direction must be plus or minus 1"))
    @boundscheck 0.0 < P.ι <= π/2.0     || throw(error("ι is outside of allowed range"))

end 



