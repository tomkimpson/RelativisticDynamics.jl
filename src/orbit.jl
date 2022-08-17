

"""
Some docstring for the run_program function
"""
function orbit(::Type{NF}=Float64;              # number format, use Float64 as default
                    kwargs...                   # all additional non-default parameters
                    ) where {NF<:AbstractFloat}


    # Setup all system parameters, universal constants etc.
    P = SystemParameters(NF=NF;kwargs...) # Parameters
    C = Constants(P)                      # Constants
    M = Model(P,C)                        # Pack all of the above into a single *Model struct 

    
    if P.model == :SphericalPhoton                       # pack all of the above into a *Model struct
        prognostic_vars = SphericalPhotonOrbit_initial_conditions(M)
    elseif P.model == :MPD
        prognostic_vars = MPD_initial_conditions(M)
    else
        println(P.model)
        println("That model selection is not defined")
        println("Please choose one of: SphericalPhoton, MPD ")
        return
    end 
    
    
        # The initial conditons. Borrowing the terminology "prognostic variables" from climate simulations
    # i.e. the xvector, pvector and svector
    #prognostic_vars = initial_conditions(M)         # initialize prognostic variables


    #Timestepping
    # Uses the solver suite DifferentialEquaitons.jl
    solution = timestepping(prognostic_vars, M)



    println("Relativistic Dynamics completed OK")
    return solution , M
end