"""
    C = Constants(P)
A struct to hold all variables which are constant over the course of the integration.
These are derived from the user-defined parameters 
"""
@with_kw struct Constants{NF<:AbstractFloat}

     
    # 1. Fundamental constants
    light_c  :: NF            # Speed of light in a vacuum, m/s
    Newton_g :: NF            # Newton's gravitational constant, m3⋅kg−1⋅s−2 
    Msolar   :: NF            # Solar mass, kg 

    # 2. Initial spacetime coordinates 
    r_initial :: NF
    θ_initial :: NF
    ϕ_initial :: NF 

    # 3. Spacetime properties 
    E :: NF  #Energy
    L :: NF  #Angular Momentum 
    Q :: NF  #Carter Constant 

    # 4. Pulsar properties
    s0 :: NF                 # Magnitude of spatial component of spin vector in natural units 
    m0 :: NF                 # Pulsar mass in natural units

    # 4. Integration properties
    Tint :: NF               # How long to integrate for

    # 5. Spin parameter
    a :: NF 

end


"""
Generator function for a Constants struct.
"""
function Constants(P::SystemParameters)

    # Fundamental constants
    light_c  = 3e8
    Newton_g = 6.67408e-11
    Msolar   = 1.989e30

    #Initial coordinates
    @unpack α = P
    r_initial = α
    θ_initial = π/2.0
    ϕ_initial = 0.0



    # Energy/Angular Momentum/Carter Constant 
    E,L,Q = ELQ(P.a,P.α,P.e,P.ι,P.orbit_dir)


    #Pulsar 
    @unpack rPSR,mPSR,mBH,p0 = P
    inertia = 0.40*mPSR*Msolar*(rPSR*1e3)^2         # Moment of inertia in SI units. Assumes a solid ball 
    convert_spin= light_c/(Newton_g*(mBH*Msolar)^2) # Multiply by this to go TO Natural units
    s0 = convert_spin*2.0*pi*inertia/p0
    m0 = mPSR/mBH

    
    #Estimate the orbital period from Kepler's 3rd. 
    #This is obviously inaccurate in relativity, but sufficient to get an approximate timescale over which to integrate 
    @unpack Norbits = P
    Tint = Norbits*2*π*α^(3/2)


    

    
   # This implies conversion to NF
    return Constants{P.NF}(light_c,Newton_g,Msolar,
                           r_initial,θ_initial,ϕ_initial,
                           E,L,Q,
                           s0,m0,
                           Tint,
                           P.a)



end






