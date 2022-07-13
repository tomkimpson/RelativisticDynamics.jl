"""Struct holding the so-called 'prognostic' variables"""
struct PrognosticVariables{NF<:AbstractFloat}
    xvector         ::AbstractVector{NF}       
end



"""Initialize a PrognosticVariables struct for an atmosphere at rest. No winds,
hence zero vorticity and divergence, but temperature, pressure and humidity are
initialised """
function initial_conditions(M::Model)

    @unpack NF,α = M.parameters




    xvector = [0.0,α,pi/2.0,0.0] # By default the starting coordinates

    

    # conversion to NF happens here
    return PrognosticVariables{NF}(xvector)
end