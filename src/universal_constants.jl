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

    # 3. Mathematical/Tensor objects
    ϵ ::AbstractArray{NF}   #levi cevita


    # 4. Integration properties
    Tint :: NF               # How long to integrate for
end


"""
Generator function for a Constants struct.
"""
function Constants(P::SystemParameters)

    println("1. Setting up universal constants")

    # Fundamental constants
    light_c  = 3e8
    Newton_g = 6.67408e-11
    Msolar   = 1.989e30

    #Initial coordinates
    @unpack α = P
    r_initial = α
    θ_initial = π/2.0
    ϕ_initial = 0.0



    # Model specific constants

    E,L,Q = ELQ(P.model,P.a,P.α,P.e,P.ι)

    println(E)
    println(L)
    println(Q)

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




    #Estimate the orbital period from Kepler's 3rd. Ofc inaccurate in relativity, but sufficient to get an approximate lengthscale
    #over which to integrate 
    @unpack Norbits = P
    Tint = Norbits*2*π*α^(3/2)

    

    
    # This implies conversion to NF
    return Constants{P.NF}(light_c,Newton_g,Msolar,
                           r_initial,θ_initial,ϕ_initial,
                           E,L,Q,
                           s0,m0,
                           levi,
                           Tint,)
end




"""
Calculate the energy, angular momentum and Carter constant for the different types of system
"""
function ELQ(model,a,α,e,ι)
    
    if model == :MPD #& e>0 

        # #ELQ circular - from Raine & Thomas 2005
        # N = sqrt(1-3/r + 2*a*r^(-3/2))
        # E = (1 - 2/r + a*r^(-3/2))/N
        # L = sqrt(r) * (1+(a/r)^2 - 2*a*r^(-3/2))/N
        # #dL = 0.1
        # #L += dL 
        # Q = 0.0

        r_periapsis= α*(1-e)
        r_apoapsis = α*(1+e)
        zm = cos(ι)

        println("peri,apo,zm")
        println(r_periapsis)
        println(r_apoapsis)
        println(zm)
        println("------------")




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
        #Need to include an ornbit direction correction here?
        D = +1
        E_numerator   = κ*ρ+2.0*η*σ-2.0*D*sqrt(σ*(σ*ϵ^2 + ρ*ϵ*κ-η*κ^2))
        E_denominator = ρ^2 + 4.0*η*σ
        E = sqrt(E_numerator/E_denominator)


        #Angular momentum 
        L = -g1*E/h1 + D*sqrt(g1^2 * E^2 +(f1*E^2 -d1)*h1)/h1



        #Carter 
        Q = zm^2*(a^2*(1.0-E^2)+L^2/(1.0-zm^2))

    else
        println(model)
        println("ELQ for that model selection is not yet defined")
        println("Please choose one of: SphericalPhoton, MPD ")
        return
    end 

    return E,L,Q


end 

# function LQ(P::SystemParameters)


#     @unpack r,a,α,e,ι,orbit_dir = P

#     if P.model == :SphericalPhoton 

#         if a == 0.0
#             println("Spherical Photon orbits do not exist for a=0. Try a non-zero value")
#             return
#         else
#             E = 1.0
#             L = - (r^3 - 3.0*r^2 + a^2*r + a^2) / (a*(r - 1.0))
#             Q = -(r^3 * (r^3 - 6.0*r^2 + 9.0 * r - 4.0*a^2))/(a^2 * (r - 1.0)^2)
#         end

#     elseif P.model == :RayTracing 


#         E = 0.0
#         L = 0.0
#         Q = 0.0

        

#         E = 0.0
#         L = 0.0
#         Q = 0.0

#     elseif P.model == :MPD


#         #ELQ circular 

#         N = 1-3/r






#         # #Some needed quantities
#         # r_periapsis= α*(1-e)
#         # println("PERIAPSIS")
#         # println(r_periapsis)

#         # zminus = cos(ι)
#         # f1 = mapping_f(r_periapsis,a,zminus)
#         # g1 = mapping_g(r_periapsis,a)
#         # h1 = mapping_h(r_periapsis,a,zminus)
#         # d1 = mapping_d(r_periapsis,a,zminus)
    
   
#         # #Angular momentum and Carter constant, taking E=1
#         # L = -g1/h1 + orbit_dir*sqrt(g1^2 + (f1-d1)*h1)/h1
#         # Q = zminus^2 * (L^2 /(1.0 - zminus^2))

    

#     else
#         println(P.model)
#         println("That model selection is not defined")
#         println("Please choose one of: SphericalPhoton, MPD ")
#         return
#     end 

#     return E,L,Q


# end 

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
