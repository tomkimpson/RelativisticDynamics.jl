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
    model::Symbol=:MPD  # Only :MPD is currently defined. 

    #BH parameters
    a::Real   = 0.1     # BH spin parameter
    mBH::Real = 4e6     # BH mass in solar masses

    #PSR parameters
    mPSR::Real = 1.4   # Pulsar mass in solar masses 
    rPSR::Real = 10.0  # Pulsar radius in km
    p0::Real   = 1e-3  # Pulsar spin period in seconds
    Sθ::Real   = pi/6  # θ angle of pulsar spin axis 
    Sϕ::Real   = 0.0   # ϕ angle of pulsar spin axis

    #Orbital parameters for the Keplerian orbit
    α::Real=50.0        # Keplerian semi major axis
    e::Real=0.10        # Keplerian eccentricity 
    ι::Real=π/6.0       # Inclination w.r.t equatorial plane in radians. The extrema of θ. 
    orbit_dir::Int=1    # Orbit direction prograde/retrograde
    Norbits = 1.0      # Number of orbits to integrate for 




end
