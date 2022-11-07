using JLD
using ChainRulesCore

function differentiate(e)


    # Load a reference solution from disk. Does not need to be differentiated
    #reference = ChainRulesCore.ignore_derivatives(load("data/example_data.jld")) 
    #println(reference)

    # Run the model. 
    println("running the model with e = ", e)
    solution,model = orbit(e=e)
    
    
    #println("the value of blob is ", blob)
    # # loss=0 # Variable must be first declared outside of ignore_derivatives to make it accesible to the scope 
    # r = 0
    # ChainRulesCore.ignore_derivatives() do # Calculate the loss w.r.t the reference



    #     p = solution(reference.t) # cast the solution onto the times of the reference
    #     r =  p[2,:]

    # end 

    # r0 = reference[2,:]
   
    # δ = r - r0
    
    # loss = sum(abs2, δ)

    r = solution[2,:]
    #println(r)
    loss = last(r) -  50.00001064382625
    #println(loss)
    #println(loss)
    #     loss = model.parameters.e#^2
    # end 

    
    #loss = model.parameters.e^2
    #loss = 0.1


    #println(loss)
    #loss = 1
    return loss

end





