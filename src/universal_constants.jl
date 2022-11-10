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
                           Tint)



end




"""
    E,L,Q = ELQ(a,α,e,ι,D)
Calculate the energy, angular momentum and Carter constant given the Keplerian orbital parameters, and the BH spin/direction
"""
function ELQ(a,α,e,ι,D)
    

    zm = cos(ι)

    r_periapsis= α*(1-e)
    r_apoapsis = α*(1+e)
    
    #Define some orbital coefficients
    f1 = mapping_f(r_periapsis,a,zm)
    g1 = mapping_g(r_periapsis,a)
    h1 = mapping_h(r_periapsis,a,zm)
    d1 = mapping_d(r_periapsis,a,zm)

    f2 = mapping_f(r_apoapsis,a,zm)
    g2 = mapping_g(r_apoapsis,a)
    h2 = mapping_h(r_apoapsis,a,zm)
    d2 = mapping_d(r_apoapsis,a,zm)

 
    #Determinants 
    κ = d1*h2 - d2*h1 
    ϵ = d1*g2 - d2*g1 
    ρ = f1*h2 - f2*h1
    η = f1*g2 - f2*g1 
    σ = g1*h2 - g2*h1 


    #Energy
    E_numerator   = κ*ρ+2.0*η*σ-2.0*D*sqrt(σ*(σ*ϵ^2 + ρ*ϵ*κ-η*κ^2))
    E_denominator = ρ^2 + 4.0*η*σ
    E = sqrt(E_numerator/E_denominator)

    #Angular momentum 
    L = -g1*E/h1 + D*sqrt(g1^2 * E^2 +(f1*E^2 -d1)*h1)/h1

    #Carter 
    Q = zm^2*(a^2*(1.0-E^2)+L^2/(1.0-zm^2))

    return E,L,Q


end 



"""
    f = mapping_f(r,a,zminus)
Mapping function `f` used when converting from Keplerian orbital parameters to constants of motion
"""
function mapping_f(r,a,zminus)
Δ = delta(r,a)
return r^4 +a^2 * (r*(r+2.0) + zminus^2 * Δ)
end


"""
    g = mapping_g(r,a)
Mapping function `g` used when converting from Keplerian orbital parameters to constants of motion
"""
function mapping_g(r,a)
return 2*a*r
end

"""
    h = mapping_h(r,a,zminus)
Mapping function `h` used when converting from Keplerian orbital parameters to constants of motion
"""
function mapping_h(r,a,zminus)
Δ = delta(r,a)
return r*(r-2.0) + (zminus^2 * Δ)/(1.0 - zminus^2)
end


"""
    d = mapping_d(r,a,zminus)
Mapping function `d` used when converting from Keplerian orbital parameters to constants of motion
"""
function mapping_d(r,a,zminus)
Δ = delta(r,a)
return (r^2 +a^2 * zminus^2)*Δ
end




