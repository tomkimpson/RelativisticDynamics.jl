

"""
Some docstring for the run_program function
"""
function orbit(::Type{NF}=Float64;         # number format, use Float64 as default
                    kwargs...                   # all additional non-default parameters
                    ) where {NF<:AbstractFloat}


    P = SystemParameters(NF=NF;kwargs...) # Parameters
    C = Constants(P)                      # Constants

    M = Model(P,C)                        # Pack all of the above into a single *Model struct 
    print(M)


    prognostic_vars = initial_conditions(M)         # initialize prognostic variables
    print(prognostic_vars)
    return "All completed OK"
end