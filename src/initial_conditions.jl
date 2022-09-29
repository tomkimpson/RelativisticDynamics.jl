"""Struct holding the so-called 'prognostic' variables"""
struct PrognosticVariables{NF<:AbstractFloat}
    xvector         ::AbstractVector{NF}       
    pvector         ::AbstractVector{NF}      
    svector         ::AbstractVector{NF}       
end




function initial_conditions(M::Model)

    # if M.parameters.model == :SphericalPhoton                       # pack all of the above into a *Model struct
    #     prognostic_vars = SphericalPhotonOrbit_initial_conditions(M)
    # elseif M.parameters.model == :RayTracing                       # pack all of the above into a *Model struct
    #         prognostic_vars = RayTracing_initial_conditions(M)
    # elseif M.parameters.model == :MPD
    #     prognostic_vars = MPD_initial_conditions(M)
    # end 
    

    if M.parameters.model == :MPD                     
        prognostic_vars = MPD_initial_conditions(M)
    else
        println("Initial conditions are not yet defined for that model")
    end 

    return prognostic_vars

end 






# """Initialize a PrognosticVariables struct for Spherical Photon orbits. 
#    Initial t,r,θ,ϕ are all set in system_parameters.jl 
#    Need to calculate pθ. see e.g. https://arxiv.org/pdf/1601.02063.pdf
#    """
# function SphericalPhotonOrbit_initial_conditions(M::Model)

#     @unpack NF,r,θ,ϕ,a = M.parameters
#     @unpack  Q,L = M.constants

#     # Initial conditions for r are set in system_parameters.jl 
#     # Initial conditions for θ, ϕ are arbitrary, also set in  system_parameters.jl 
#     # We now calculate the initial value of pθ
    
#     Σ = sigma(r,θ,a)
#     θdot2 = (Q - (L^2 / sin(θ)^2 - a^2)*cos(θ)^2)/Σ^2
#     θdot = sqrt(θdot2) # plus or minus?
#     pθ = Σ*θdot

#     # 4- position

#     xvector = [0.0,r,θ,ϕ]     # By default the starting coordinates
#     pvector = [0.0,0.0,pθ,0.0]
#     svector = [0.0,0.0,0.0,0.0] # dont track spin for these guys


#     # conversion to NF happens here
#     return PrognosticVariables{NF}(xvector,pvector,svector)
# end


# """Initialize a PrognosticVariables struct for Spherical Photon orbits. 
#    Initial t,r,θ,ϕ are all set in system_parameters.jl 
#    Need to calculate pθ. see e.g. https://arxiv.org/pdf/1601.02063.pdf
#    """
# function RayTracing_initial_conditions(M::Model)

#     @unpack NF,r,θ,ϕ,a = M.parameters
#     @unpack  Q,L = M.constants

#     println("CALLED RAY TRACING INITIaL CONDITONS")



#     α_c = 7
#     β_c = 0

#     θ_obs = π/2.0
#     r_obs = 1000

#     xp = sqrt(r_obs^2 + a^2)*sin(θ_obs) - β_c*cos(θ_obs)
#     yp = α_c 
#     zp = r_obs*cos(θ_obs) + β_c*sin(θ_obs)

#     #Convert to BL coords 
#     w = xp^2 + yp^2 + zp^2 -a^2
#     rbar = sqrt((w+sqrt(w^2 + 4*a^2*zp))/2)
#     θbar = acos(zp/rbar)
#     ϕbar = atan2(yp,xp)

#     Σ = sigma(rbar,θbar,a)
#     Δ = delta(rbar,a)

#     u = sqrt(rbar^2 +a^2)
#     v = -sin(θ_obs)*cos(ϕbar)
#     zdot = -1.0

#     rdot = -zdot*(-u^2*cos(θ_obs)*cos(θbar)+rbar*u*v*sin(θbar))/Σ


#     beta_coord = 0.0



    
#     Σ = sigma(r,θ,a)
#     θdot2 = (Q - (L^2 / sin(θ)^2 - a^2)*cos(θ)^2)/Σ^2
#     θdot = sqrt(θdot2) # plus or minus?
#     pθ = Σ*θdot

#     # 4- position

#     xvector = [0.0,r,θ,ϕ]     # By default the starting coordinates
#     pvector = [0.0,0.0,pθ,0.0]
#     svector = [0.0,0.0,0.0,0.0] # dont track spin for these guys


#     # conversion to NF happens here
#     return PrognosticVariables{NF}(xvector,pvector,svector)
# end






"""Setup the initial conditions for the MPD orbital dynamics"""
function MPD_initial_conditions(M::Model)

    println("Setting up initial conditions")

   @unpack NF = M.parameters
   # @unpack L,Q,s0,m0 = M.constants 



    # 1. Four- position
    @unpack r_initial,θ_initial, ϕ_initial = M.constants
    xvector = [0.0,r_initial,θ_initial, ϕ_initial] # Starting coordinates



    #1.1 Define some useful quantities
    @unpack a = M.parameters
    r,θ= xvector[2],xvector[3] 
    println("Initial coordinates are:")
    println(r)
    println(θ)

    println("spin")
    println(a)

    Δ = delta(r,a)
    Σ = sigma(r,θ,a)
    g  = covariant_metric(xvector,a)

   
    # 2. Four - momentum 
    @unpack E,L,Q,m0 = M.constants 


    #NEED TO UPDATE THESE FOR E
    Pbar = E*(r^2+a^2) - a*L 
    Tbar = (r^2+a^2)*Pbar/Δ -a*(a*E*sin(θ)^2-L)
    Rbar = ((r^2+a^2)*E  -a*L)^2 -Δ*(r^2 + (L-a*E)^2+Q)
    θbar = Q - ((1-E^2)*a^2 + L^2/sin(θ)^2)*cos(θ)^2
    ϕbar = a*Pbar/Δ -a*E + L/sin(θ)^2

    println("RBAR=")
    println(Rbar)

    println("thet bar")
    println(θbar)

    println(θ)
    println(sin(θ))

    #These are 4 velocities from Schmidt 2002.
    #Initial Rdot is +ve as standard
    tdot = Tbar/Σ
    rdot = sqrt(Rbar)/Σ
    θdot = sqrt(θbar)/Σ
    ϕdot = ϕbar/Σ

   

    #Turn 4 velocity into 4 momentum
    pvector = m0*[tdot,rdot,θdot,ϕdot]
    pvector_covar = convert_to_covariant(g,pvector)


    #3. Four - spin
    @unpack Sθ,Sϕ = M.parameters
    @unpack s0 = M.constants 

    svector = [0.0,0.0,0.0,0.0]

    #Set the spatial components of the spin vector
    svector[2] = s0 * sin(Sθ) * cos(Sϕ)
    svector[3] = -s0 *cos(Sθ)/r 
    svector[4] = s0*sin(Sθ)*sin(Sϕ)/(r*sin(θ)) 
    svector[1] = -(svector[2]*pvector_covar[2] + svector[3]*pvector_covar[3]+svector[4]*pvector_covar[4])/pvector_covar[1] # Enforces the spin condition


    # conversion to NF happens here
    return PrognosticVariables{NF}(xvector,pvector,svector)
end