

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

    
    prognostic_vars = initial_conditions(M)
    solution = timestepping(prognostic_vars, M)
    

    
    


    #Timestepping
    # Uses the solver suite DifferentialEquaitons.jl
    #println("Starting the timestepping")
    #solution = timestepping(prognostic_vars, M)
    #println("Relativistic Dynamics completed OK")
    #solution = 1
    #display(solution)
    return solution , M
end