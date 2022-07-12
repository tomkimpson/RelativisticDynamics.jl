"""
    C = Constants(P)
A struct to hold all variables which are constant over the course of the integration.
"""
@with_kw struct Constants{NF<:AbstractFloat}

    # Fundamental constants 
    light_c :: NF

end


"""
Generator function for a Constants struct.
"""
function Constants(P::SystemParameters)

    # PHYSICAL CONSTANTS
    @unpack a = P

    light_c = a+1.0
    

    # This implies conversion to NF
    return Constants{P.NF}(light_c)
end