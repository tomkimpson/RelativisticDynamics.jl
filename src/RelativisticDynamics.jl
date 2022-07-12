module RelativisticDynamics




# Imports
import Parameters: @with_kw, @unpack
#import ProtoStructs: @proto





# Exports
export  orbit, SystemParameters, Constants

#Includes
include("system_parameters.jl")
include("universal_constants.jl")
include("orbit.jl")





print("Welcome to the Relativistic Dynamics module")
end
