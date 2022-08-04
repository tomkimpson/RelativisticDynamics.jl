"""
    C = Constants(P)
A struct to hold all variables which are constant over the course of the integration.
"""
@with_kw struct Constants{NF<:AbstractFloat}

     
    # Fundamental constants
    light_c  :: NF            # Speed of light in a vacuum, m/s
    Newton_g :: NF            # Newton's gravitational constant, m3⋅kg−1⋅s−2 
    Msolar   :: NF            # Solar mass, kg 


    # Constants derived from user defined parameters 
    E   :: NF    
    L   :: NF   
    Q   :: NF      
    
    
    Φ :: NF 
    u0:: NF
    u1 ::NF 


    # Pulsar constants
    s0 :: NF # Magnitude of spatial component of spin vector in natural units 


    # BH constants
    rH :: NF #Black hole outer horizon radius

end


"""
Generator function for a Constants struct.
"""
function Constants(P::SystemParameters)

    # Unpack system parameters
    @unpack a,ι,α,e,orbit_dir,mPSR,rPSR,p0,mBH,r = P

    # Fundamental constants
    light_c  = 3e8
    Newton_g = 6.67408e-11
    Msolar   = 1.989e30


    #SphericalPhoton orbit constants 
    # Define constants of motion 
    #Σ = sigma(r,θ,a)
    Φ = - (r^3 - 3.0*r^2 + a^2*r + a^2) / (a*(r - 1.0))
    Q = -(r^3 * (r^3 - 6.0*r^2 + 9.0 * r - 4.0*a^2))/(a^2 * (r - 1.0)^2)


    # Define u0 function 
    u0 = ((a^2 - Q - Φ^2) + sqrt((a^2 - Q - Φ^2)^2 + 4.0*a^2 * Q))/(2.0*a^2)
    u1 = ((a^2 - Q - Φ^2) - sqrt((a^2 - Q - Φ^2)^2 + 4.0*a^2 * Q))/(2.0*a^2)


    #Constants/Orbit mapping 
    zminus = sin(ι)^2


    r_periapsis= α*(1-e)
    r_apoapsis = α*(1+e)

    f1 = mapping_f(r_periapsis,a,zminus)
    g1 = mapping_g(r_periapsis,a)
    h1 = mapping_h(r_periapsis,a,zminus)
    d1 = mapping_d(r_periapsis,a,zminus)

    f2 = mapping_f(r_apoapsis,a,zminus)
    g2 = mapping_g(r_apoapsis,a)
    h2 = mapping_h(r_apoapsis,a,zminus)
    d2 = mapping_d(r_apoapsis,a,zminus)


    κ = d1*h2 - d2*h1
    ϵ = d1*g2 - d2*g1
    ρ = f1*h2 - f2*h1
    η = f1*g2 - f2*g1
    σ = g1*h2 - g2*h1


    #Energy
    E2 = κ*ρ + 2.0*ϵ*σ + orbit_dir*2.0*sqrt(σ*(σ*ϵ^2 + ρ*ϵ*κ - η*κ^2)) / (ρ^2 + 4.0*η*ρ)
    E = sqrt(E2)
 
    #Angular Momentum 
    L = -g1*E/h1 + orbit_dir*sqrt(g1^2*E2 + (f1*E2-d1)*h1)/h1


    #Carter Constant 
    #Q = zminus * (a^2 * (1.0 - E2 ) + L^2 /(1.0 - zminus))

    


    #Pulsar 
    inertia = 0.40*mPSR*Msolar*(rPSR*1e3)^2 # Moment of inertia in SI units. Assumes solid ball 
    convert_spin= light_c/(Newton_g*(mBH*Msolar)^2) # Multiply by this to go TO Natural units
    s0 = convert_spin*2.0*pi*inertia/p0


    #Black hole 
    rH = 1.0 + sqrt(1.0 - a^2)




    # This implies conversion to NF
    return Constants{P.NF}(light_c,Newton_g,Msolar,
                           E,L,Q,
                           Φ,u0,u1,
                           s0,
                           rH)
end





function EQL(P::SystemParameters)



end 