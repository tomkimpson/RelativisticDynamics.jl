

function differentiate(::Type{NF}=Float64;              # number format, use Float64 as default
                    kwargs...                   # all additional non-default parameters
                    ) where {NF<:AbstractFloat}

    reference_solution,reference_model = orbit() #Always uses the defaults

    solution,model = orbit(NF=NF;kwargs...)      #uses the user-passed arguments

  
end






# function loss_function(::Type{NF}=Float64;              # number format, use Float64 as default
#                        kwargs...                   # all additional non-default parameters
#                       ) where {NF<:AbstractFloat}


# true_solution,true_model = orbit()


# new_solution,new_model = orbit(NF=NF;kwargs...)


# end 