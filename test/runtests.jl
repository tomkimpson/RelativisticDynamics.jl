using RelativisticDynamics
using Test,Documenter
using LinearAlgebra
using Distributions
using Tullio


doctest(RelativisticDynamics)

# GENERAL
include("metric.jl")
include("useful_functions.jl")
include("initial_conditions.jl")
include("universal_constants.jl")
include("initial_conditions.jl")
include("timestepping.jl")
include("number_formats.jl")
include("orbit.jl")