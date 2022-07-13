

"""
Some docstring for the run_program function
"""
function orbit(::Type{NF}=Float64;         # number format, use Float64 as default
                    kwargs...                   # all additional non-default parameters
                    ) where {NF<:AbstractFloat}


    P = SystemParameters(NF=NF;kwargs...) # Parameters
    C = Constants(P)                      # Constants
    #metric = initialize_metric(P)         # Initialize the array to hold the metric 
    M = Model(P,C)                 # Pack all of the above into a single *Model struct 


    prognostic_vars = initial_conditions(M)         # initialize prognostic variables
    return "All completed OK"
end