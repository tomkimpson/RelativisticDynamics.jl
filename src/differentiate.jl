using JLD
using ChainRulesCore
function differentiate(reference, ::Type{NF}=Float64;              # number format, use Float64 as default
                        kwargs...                   # all additional non-default parameters
                        ) where {NF<:AbstractFloat}



    # Load a reference solution from disk. Does not need to be differentiated
    #reference = ChainRulesCore.ignore_derivatives(load("data/example_data.jld")) 
    #println(reference)

    # Run the model. 
    solution,model = orbit(NF=NF;kwargs...)
    
    
    # # loss=0 # Variable must be first declared outside of ignore_derivatives to make it accesible to the scope 
    # r = 0
    # ChainRulesCore.ignore_derivatives() do # Calculate the loss w.r.t the reference



    #     p = solution(reference.t) # cast the solution onto the times of the reference
    #     r =  p[2,:]

    # end 

    # r0 = reference[2,:]
   
    # δ = r - r0
    
    # loss = sum(abs2, δ)

    r = solution[2,33]
    println(r)
    #loss = sum(r)
    #println(loss)
    #println(loss)
    #     loss = model.parameters.e#^2
    # end 

    
    #loss = model.parameters.e^2
    #loss = 0.1


    #println(loss)
    loss = r
    return loss

end






# function loss_function(::Type{NF}=Float64;              # number format, use Float64 as default
#                        kwargs...                   # all additional non-default parameters
#                       ) where {NF<:AbstractFloat}


# true_solution,true_model = orbit()


# new_solution,new_model = orbit(NF=NF;kwargs...)


# end 