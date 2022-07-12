"""
    P = Parameters(kwargs...)
A struct to hold all model parameters that may be changed by the user.
The struct uses keywords such that default values can be changed at creation.
The default values of the keywords define the default model setup.
"""
@with_kw struct SystemParameters

    # NUMBER FORMATS
    NF::DataType       # Number format. Default is defined in orbit.jl

    #BH parameters
    a::Real   = 0.9     # BH spin parameter
    mBH::Real = 4e6     # BH mass in solar masses

    #PSR parameters
    mPSR::Real = 1.4   # Pulsar mass in solar masses 
    rPSR::Real = 10.0  # Pulsar radius in km
    p0::Real   = 1e-3  # Pulsar spin period in seconds
    Sθ::Real   = pi/6  # θ angle of pulsar spin axis 
    Sϕ::Real   = 0.0   # \phi angle of pulsar spin axis

    #Orbital parameters
    α::Real=30.0       # Keplerian semi major axis

end


#Could also use protostructs whilst in dev: https://github.com/BeastyBlacksmith/ProtoStructs.jl
# @proto struct SystemParameters
#     NF::DataType       # Number format. Default is defined in orbit.jl
#     a::Int = 1
#     b::Float64 = 2.0
# end
