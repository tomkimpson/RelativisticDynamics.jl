abstract type ModelSetup end

"""
"""
struct Model{NF<:AbstractFloat} <: ModelSetup
    parameters::SystemParameters
    constants::Constants{NF}
end

