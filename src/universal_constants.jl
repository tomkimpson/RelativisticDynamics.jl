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
    L   :: NF                # Angular Momentum 
    Q   :: NF                # Carter Constant 
    
    # 2.2 Pulsar properties
    s0 :: NF                 # Magnitude of spatial component of spin vector in natural units 
    m0 :: NF                 # Pulsar mass in natural units

    # 3. Mathematical/Tensor objects
    ϵ ::AbstractArray{NF}   #levi cevita

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
    L,Q = LQ(P)
    
    

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
                    levi[i,j,k,l] = levicivita(permutation_vector) #This is [i,j,k,l] from e.g. https://mathworld.wolfram.com/PermutationTensor.html
                end
            end 
        end
    end 

    
    # This implies conversion to NF
    return Constants{P.NF}(light_c,Newton_g,Msolar,
                           L,Q,
                           s0,m0,
                           levi)
end




"""
Calculate the energy, angular momentum and Carter constant for the different types of system
"""
function LQ(P::SystemParameters)


    @unpack r,a,α,e,ι,orbit_dir = P

    if P.model == :SphericalPhoton 

        if a == 0.0
            println("Spherical Photon orbits do not exist for a=0. Try a non-zero value")
            return
        else
            L = - (r^3 - 3.0*r^2 + a^2*r + a^2) / (a*(r - 1.0))
            Q = -(r^3 * (r^3 - 6.0*r^2 + 9.0 * r - 4.0*a^2))/(a^2 * (r - 1.0)^2)
        end
    elseif P.model == :MPD

        #Some needed quantities
        r_periapsis= α*(1-e)
        zminus = cos(ι)
        f1 = mapping_f(r_periapsis,a,zminus)
        g1 = mapping_g(r_periapsis,a)
        h1 = mapping_h(r_periapsis,a,zminus)
        d1 = mapping_d(r_periapsis,a,zminus)
    
   
        #Angular momentum and Carter constant, taking E=1
        L = -g1/h1 + orbit_dir*sqrt(g1^2 + (f1-d1)*h1)/h1
        Q = zminus^2 * (L^2 /(1.0 - zminus^2))

    

    else
        println(P.model)
        println("That model selection is not defined")
        println("Please choose one of: SphericalPhoton, MPD ")
        return
    end 

    return L,Q


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
