"""Struct holding the so-called 'prognostic' variables"""
struct PrognosticVariables{NF<:AbstractFloat}
    xvector         ::AbstractVector{NF}       
    pvector         ::AbstractVector{NF}      
    svector         ::AbstractVector{NF}       
end


"""Initialize a PrognosticVariables struct for an atmosphere at rest. No winds,
hence zero vorticity and divergence, but temperature, pressure and humidity are
initialised """
function SphericalPhotonOrbit_initial_conditions(M::Model)

    @unpack NF,r,θ,ϕ,a = M.parameters
    @unpack  Q,L = M.constants

    # Initial conditions for r are set in system_parameters.jl 
    # Initial conditions for θ, ϕ are arbitrary, also set in  system_parameters.jl 
    # We now calculate the initial value of pθ
    
    Σ = sigma(r,θ,a)
    θdot2 = (Q - (L^2 / sin(θ)^2 - a^2)*cos(θ)^2)/Σ^2
    θdot = sqrt(θdot2) # plus or minus?
    pθ = Σ*θdot

    # 4- position
    #
    #χ = acos(cos(θ)/sqrt(u0))
    xvector = [0.0,r,θ,ϕ]     # By default the starting coordinates
    pvector = [0.0,0.0,pθ,0.0]
    svector = [0.0,0.0,0.0,0.0] #dont track spin or momentum for these guys


    # conversion to NF happens here
    return PrognosticVariables{NF}(xvector,pvector,svector)
end

"""Initialize a PrognosticVariables struct for an atmosphere at rest. No winds,
hence zero vorticity and divergence, but temperature, pressure and humidity are
initialised """
function MPD_initial_conditions(M::Model)

    @unpack NF,α,θ,ϕ,a,mPSR,mBH,Sθ,Sϕ = M.parameters
    @unpack E,L,Q,s0,m0 = M.constants 

    println("Welcome to MPD initial conditons")

    # 4- position
    xvector = [0.0,α,θ,ϕ] # By default the starting coordinates


    r,θ = xvector[2],xvector[3] #Get the initial r and theta 

    Δ = delta(r,a)
    Σ = sigma(r,θ,a)


    # Metric 
    g  = covariant_metric(xvector,a)
    metric_contra = contravariant_metric(xvector,a)


    println([E,L,Q])

    # 4 - momentum 
    PP = E*(r^2+a^2) - a*L 


    TT = (r^2+a^2)*PP/Δ -a*(a*E*sin(θ)^2-L)
    RR = ((r^2+a^2)*E  -a*L)^2 -Δ*(r^2 + (L-a*E)^2+Q)
    ThTh = Q - ((1.0-E^2)*a^2 + L^2/sin(θ)^2)*cos(θ)^2
    PhPh = a*PP/Δ -a*E + L/sin(θ)^2



    # T = (r^2 + a^2)*(E*(r^2 + a^2) -a*L)/Δ -a*(a*E*sin(θ)^2 - L)
    # R = ((r^2 + a^2)*E -a*L)^2 - Δ*(r^2 + (L -a*E) + Q)
    # Θ = Q - ((1.0 - E^2)*a^2 + L^2/sin(θ)^2)*cos(θ)^2                                                 # Note that this is a capital \Theta, not a lower case \theta
    # Φ = a*(E*(r^2 + a^2) -a*L)/Δ -a*E + L/sin(θ)^2 



    # tdot = T/Σ
    # rdot = -sqrt(R)/Σ
    # θdot = -sqrt(Θ)/Σ
    # ϕdot = Φ/Σ




    #pvector = m0*[tdot,rdot,θdot,ϕdot]





    tdot = -TT/Σ
    rdot = -sqrt(RR)/Σ
    θdot = -sqrt(ThTh)/Σ
    ϕdot = -PhPh/Σ



    uvector = [tdot,rdot,θdot,ϕdot]



    # @tensor begin
    #     dx[α] := -pvector[α]#(pvector[α]+0.50*Stensor[α,β]*Riemann[β,γ,μ,ν]*pvector[γ]*Stensor[μ,ν]/(m0^2 + division_scalar))/m0^2
    # end

    @tensor begin 
        Vsq = g[μ,ν]*uvector[μ]*uvector[ν]
    end 


    println("VSQ")
    println(Vsq)
    PV = -sqrt(-1.0/Vsq)






    pvector_covar = convert_to_covariant(g,pvector)


    # 4 - spin
    svector = [0.0,0.0,0.0,0.0]

    #Set the spatial components of the spin vector
    svector[2] = s0 * sin(Sθ) * cos(Sϕ)
    svector[3] = -s0 *cos(Sθ)/r 
    svector[4] = s0*sin(Sθ)*sin(Sϕ)/(r*sin(θ)) 
    svector[1] = -(svector[2]*pvector_covar[2] + svector[3]*pvector_covar[3]+svector[4]*pvector_covar[4])/pvector_covar[1]


    # conversion to NF happens here
    return PrognosticVariables{NF}(xvector,pvector,svector)
end