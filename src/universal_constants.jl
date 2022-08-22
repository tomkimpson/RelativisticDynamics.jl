"""
    C = Constants(P)
A struct to hold all variables which are constant over the course of the integration.
"""
@with_kw struct Constants{NF<:AbstractFloat}

     
    # 1. Fundamental constants
    light_c  :: NF            # Speed of light in a vacuum, m/s
    Newton_g :: NF            # Newton's gravitational constant, m3⋅kg−1⋅s−2 
    Msolar   :: NF            # Solar mass, kg 


    # 2. Constants derived from user defined parameters 
    # 2.1 Spacetime properties
    E   :: NF                # Energy
    L   :: NF                # Angular Momentum 
    Q   :: NF                # Carter Constant 
    
    # 2.2 Pulsar properties
    s0 :: NF                 # Magnitude of spatial component of spin vector in natural units 
    m0 :: NF                 # Pulsar mass in natural units

    # 3. Mathematical/Tensor objects
    ϵ ::AbstractArray{NF}   #levi cevit

end


"""
Generator function for a Constants struct.
"""
function Constants(P::SystemParameters)

    # Unpack system parameters
    #@unpack a,ι,α,e,orbit_dir,mPSR,rPSR,p0,mBH,r = P

    # Fundamental constants
    light_c  = 3e8
    Newton_g = 6.67408e-11
    Msolar   = 1.989e30



    # Model specific constants
    E,L,Q = ELQ(P)
    
    

    #Pulsar 
    @unpack rPSR,mPSR,mBH,p0 = P
    inertia = 0.40*mPSR*Msolar*(rPSR*1e3)^2         # Moment of inertia in SI units. Assumes solid ball 
    convert_spin= light_c/(Newton_g*(mBH*Msolar)^2) # Multiply by this to go TO Natural units
    s0 = convert_spin*2.0*pi*inertia/p0
    m0 = mPSR/mBH


    # Levi civita tensor
    levi = zeros(Float64,4,4,4,4) 
    for i in 1:4
        for j in 1:4
            for k in 1:4
                for l in 1:4
                    permutation_vector = [i,j,k,l]
                    levi[i,j,k,l] = levicivita(permutation_vector)
                end
            end 
        end
    end 

    
    # This implies conversion to NF
    return Constants{P.NF}(light_c,Newton_g,Msolar,
                           E,L,Q,
                           s0,m0,
                           levi)
end




"""
Calculate the energy, angular momentum and Carter constant for the different types of system
"""
function ELQ(P::SystemParameters)


    #@unpack a,ι,α,e,orbit_dir,mPSR,rPSR,p0,mBH,r = P

    @unpack r,a,α,e,ι,orbit_dir = P


    if P.model == :SphericalPhoton 

        E = 1.0
        L = - (r^3 - 3.0*r^2 + a^2*r + a^2) / (a*(r - 1.0))
        Q = -(r^3 * (r^3 - 6.0*r^2 + 9.0 * r - 4.0*a^2))/(a^2 * (r - 1.0)^2)

    elseif P.model == :MPD

        r_periapsis= α*(1-e)
        r_apoapsis = α*(1+e)
        zminus = cos(ι)
    
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
        E2 = κ*ρ + 2.0*ϵ*σ + orbit_dir*2.0*sqrt(σ*(σ*ϵ^2 + ρ*ϵ*κ - η*κ^2)) / (ρ^2 + 4.0*η*σ)
        E = sqrt(E2)
     
        #Angular Momentum 
        L = -g1*E/h1 + orbit_dir*sqrt(g1^2*E2 + (f1*E2-d1)*h1)/h1
        #L = (ρ*E2 - κ)/(2.0*E*σ)
    
        #Carter Constant 
        Q = zminus^2 * (a^2 * (1.0 - E2 ) + L^2 /(1.0 - zminus^2))
    
        #Keplerian expressions for reference. See Schmidt 2002 appendix 
        # SLR = α*(1-e^2)
        # zm = cos(ι)
        # E = sqrt(1.0 - (1-e^2)/SLR)
        # L = sqrt((1-zm^2)*SLR)
        # Q = zm^2*SLR

        # Useful numbers
        # E = 0.999390486044721
        # L = 14.515059110545808
        # Q = 210.68716034079137

    else
        println(P.model)
        println("That model selection is not defined")
        println("Please choose one of: SphericalPhoton, MPD ")
        return
    end 

    return E,L,Q


end 