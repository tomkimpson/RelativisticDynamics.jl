module RelativisticDynamics




# Imports
import Parameters: @with_kw, @unpack
import DifferentialEquations
using Tullio,Combinatorics,LinearAlgebra
using Zygote
# Exports
export  orbit, PlotTrajectory, StackedPlot

#Includes
include("system_parameters.jl")
include("useful_functions.jl")
include("metric.jl")
include("universal_constants.jl")
include("initial_conditions.jl")
include("timestepping.jl")
include("orbit.jl")
include("plotting.jl")


print("Welcome to the Relativistic Dynamics module!")
end
