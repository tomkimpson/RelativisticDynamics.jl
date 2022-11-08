

function differentiate(reference,e)

        # Run the model using the new value of e
        solution,model = orbit(e=e)
        
    
        r0 = last(reference[2,:])  # The last r value of the reference solution
        r = last(solution[2,:])    # The last r value of the new solution
        loss = abs(r0 - r)         # Loss
      
        return loss
    

end 