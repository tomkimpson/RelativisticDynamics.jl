

"""
Some docstring for the run_program function"""
function orbit(::Type{NF}=Float64;         # number format, use Float64 as default
                    kwargs...                   # all additional non-default parameters
                    ) where {NF<:AbstractFloat}


    P = SystemParameters(NF=NF;kwargs...) # Parameters
    C = Constants(P)                      # Constants
    print(P)
    print(C)
    return "All completed OK"
end