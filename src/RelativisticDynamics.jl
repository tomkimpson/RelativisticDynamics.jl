module RelativisticDynamics




# Imports
import Parameters: @with_kw, @unpack
#import ProtoStructs: @proto





# Exports
export  orbit, SystemParameters, Constants, PrognosticVariables

#Includes
include("system_parameters.jl")
include("useful_functions.jl")
include("universal_constants.jl")
include("model.jl")
include("initial_conditions.jl")
include("orbit.jl")





print("Welcome to the Relativistic Dynamics module")
end
