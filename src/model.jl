abstract type ModelSetup end

"""
    M = Model(P,C) 
The model struct which holds all the parameters (P) and constants (C)
"""
struct Model{NF<:AbstractFloat} <: ModelSetup
    parameters::SystemParameters
    constants::Constants{NF}
end

