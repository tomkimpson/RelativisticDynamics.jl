module RelativisticDynamics




# Imports
import Parameters: @with_kw, @unpack

import DifferentialEquations
using TensorOperations, Einsum #, Symbolics



# Exports
export  orbit, SystemParameters, Constants, PrognosticVariables,PlotTrajectory,AnimateTrajectory

#Includes
include("system_parameters.jl")
include("useful_functions.jl")
include("metric.jl")
include("universal_constants.jl")
include("model.jl")
include("initial_conditions.jl")
include("timestepping.jl")
include("orbit.jl")





print("Welcome to the Relativistic Dynamics module!")
end
