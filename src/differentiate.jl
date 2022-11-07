using JLD
using ChainRulesCore
using Zygote

function differentiate(reference,e)


    # Load a reference solution from disk. Does not need to be differentiated
    #reference = ChainRulesCore.ignore_derivatives(load("data/example_data.jld")) 
    #println(reference)

    # Run the model. 
    println("running the model with e = ", e)
    solution,model = orbit(e=e)
    
    r0 = last(reference[2,:])
    r = last(solution[2,:])
    loss = abs(r0 - r)
  

    return loss

end







function difference()

    reference_solution,reference_model = orbit()
    println("---------------------------------")

    N = 5
    data = zeros(Float64,N,2)
    for (i,v) in pairs(range(0.1,0.9,N))

        loss = Zygote.gradient(x -> differentiate(reference_solution,x),v)[1]
        #loss = differentiate(reference_solution,v)
        data[i,1] = v 
        data[i,2] = loss

    end


    return data
end 












#SCRACTH




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





    # r = solution[2,:]
    # #println(r)
    # #loss = last(r) -  48.91429371658843
    # loss = model.constants.L - 3.6283017608788204
    # #loss = model.parameters.e - 0.10
    # println("loss = ",)
    # println(loss)
    # println(last(r))
    # #println(loss)
    # #loss = model.parameters.e^2
    # # end 

    
    # #loss = model.parameters.e^2
    # #loss = 0.1


    # #println(loss)
    # #loss = 1
    # #println(loss)