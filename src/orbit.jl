

"""
The main output function from this package.
"""
function orbit(::Type{NF}=Float64;              # number format, use Float64 as default
                    kwargs...                   # all additional non-default parameters
                    ) where {NF<:AbstractFloat}


    # Setup all system parameters, universal constants etc.
    P = SystemParameters(NF=NF;kwargs...) # Parameters
    bounds_checks(P)

    C = Constants(P)                      # Constants
    M = Model(P,C)                        # Pack all of the above into a single *Model struct 

    #Initial conditions 
    initialization = initial_conditions(M)

    #Evolve in time
    solution = timestepping(initialization, M)
    
    println("Completed OK")
    return solution , M
end



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