"""
    P = Parameters(kwargs...)
A struct to hold all model parameters that may be changed by the user.
The struct uses keywords such that default values can be changed at creation.
The default values of the keywords define the default model setup.
"""
@with_kw struct SystemParameters

    # NUMBER FORMATS
    NF::DataType       # Number format. Default is defined in orbit.jl


    #Type of system to integrate 
    model::Symbol=:SphericalPhoton           # :barotropic, :shallowwater, or :primitive

    #BH parameters
    a::Real   = 0.9     # BH spin parameter
    mBH::Real = 4e6     # BH mass in solar masses

    #PSR parameters
    mPSR::Real = 1.4   # Pulsar mass in solar masses 
    rPSR::Real = 10.0  # Pulsar radius in km
    p0::Real   = 1e-3  # Pulsar spin period in seconds
    Sθ::Real   = pi/6  # θ angle of pulsar spin axis 
    Sϕ::Real   = 0.0   # \phi angle of pulsar spin axis

    #Orbital parameters for the Keplerian orbit
    α::Real=300.0       # Keplerian semi major axis
    e::Real=0.1         # Keplerian eccentricity 
    ι::Real=pi/12       # Inclination w.r.t equatorial plane in radians 
    orbit_dir::Int=1    # Orbit direction prograde/retrograde

    #Orbital parameters for the Spherical Photon orbit 
    rmin = 2.0 * (1.0 + cos(2.0/3.0 * acos(-a)))
    rmax = 2.0 * (1.0 + cos(2.0/3.0 * acos(a)))
    r::Real=(rmin+rmax)/2.0 
    θ::Real = π/6.0
    ϕ::Real = 0.0
    Tint::Real = 100.0 # How long to integrate for 


end


#Could also use protostructs whilst in dev: https://github.com/BeastyBlacksmith/ProtoStructs.jl
# @proto struct SystemParameters
#     NF::DataType       # Number format. Default is defined in orbit.jl
#     a::Int = 1
#     b::Float64 = 2.0
# end
