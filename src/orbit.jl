
"""
    solution,model = orbit(NF,kwargs...)
    
Runs RelativisticDynamics.jl with number format `NF` and any additional parameters in the keyword arguments
`kwargs...`. Any unspecified parameters will use the default values as defined in `src/system_parameters.jl`."""
function orbit(::Type{NF}=Float64;              # number format, use Float64 as default
               kwargs...                        # all additional non-default parameters
               ) where {NF<:AbstractFloat}

    # Setup all system parameters, universal constants etc.
    P = SystemParameters(NF=NF;kwargs...) # Parameters
    bounds_checks(P)                      # Check all parameters are reasonable
    C = Constants(P)                      # Constants
    M = Model(P,C)                        # Pack all of the above into a single *Model struct 

    #Initial conditions 
    initialization = initial_conditions(M)

    #Evolve in time
    solution = timestepping(initialization, M)
    return solution, M

end


"""
    bounds_checks(P)
Look before you leap - are the user-specified kwargs in parameters physical and reasonable?
"""
function bounds_checks(P::SystemParameters)
    @boundscheck 0.0<=P.a<1.0           || throw(error("Spin parameter a is out of bounds"))
    @boundscheck P.mBH/P.mPSR >= 1e3    || throw(error("Mass ratio is too small"))
    @boundscheck 1 <= P.rPSR <= 100     || throw(error("Pulsar radius is unphysical"))
    @boundscheck 0.0< P.e<1.0          || throw(error("Eccentricity is outside range"))
    @boundscheck P.orbit_dir in [-1,1]  || throw(error("Orbit direction must be plus or minus 1"))
    @boundscheck 0.0 < P.ι <= π/2.0     || throw(error("ι is outside of allowed range"))

end 


"""
    M = Model(P,C) 
The model struct which holds all the parameters (P) and constants (C)
"""
struct Model{NF<:AbstractFloat} #<: ModelSetup
    parameters::SystemParameters
    constants::Constants{NF}
end

