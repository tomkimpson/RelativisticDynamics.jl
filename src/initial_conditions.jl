"""Struct holding the so-called 'prognostic' variables"""
struct PrognosticVariables{NF<:AbstractFloat}
    xvector         ::AbstractVector{NF}       
    pvector         ::AbstractVector{NF}      
    svector         ::AbstractVector{NF}       
end




function initial_conditions(M::Model)

    if M.parameters.model == :MPD                     
        prognostic_vars = MPD_initial_conditions(M)
    else
        println("Initial conditions are not yet defined for that model")
    end 

    return prognostic_vars

end 





"""Setup the initial conditions for the MPD orbital dynamics"""
function MPD_initial_conditions(M::Model)


    @unpack NF = M.parameters
   
    # 1. Four- position
    @unpack r_initial,θ_initial, ϕ_initial = M.constants
    xvector = [0.0,r_initial,θ_initial, ϕ_initial] # Starting coordinates



    #1.1 Define some useful quantities
    @unpack a = M.parameters
    r,θ= xvector[2],xvector[3] 
    Δ = delta(r,a)
    Σ = sigma(r,θ,a)
    g  = covariant_metric(xvector,a)

   
    # 2. Four - momentum 
    @unpack E,L,Q,m0 = M.constants 


    #These are 4 velocities from Schmidt 2002.
    #Initial Rdot is +ve as standard
    Pbar = E*(r^2+a^2) - a*L 
    Tbar = (r^2+a^2)*Pbar/Δ -a*(a*E*sin(θ)^2-L)
    Rbar = ((r^2+a^2)*E  -a*L)^2 -Δ*(r^2 + (L-a*E)^2+Q)
    θbar = Q - ((1-E^2)*a^2 + L^2/sin(θ)^2)*cos(θ)^2
    ϕbar = a*Pbar/Δ -a*E + L/sin(θ)^2

    
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
    svector[2] =  s0 * sin(Sθ) * cos(Sϕ)
    svector[3] = -s0 *cos(Sθ)/r 
    svector[4] =  s0*sin(Sθ)*sin(Sϕ)/(r*sin(θ)) 
    svector[1] = -(svector[2]*pvector_covar[2] + svector[3]*pvector_covar[3]+svector[4]*pvector_covar[4])/pvector_covar[1] # Enforces the spin condition


    # conversion to NF happens here
    return PrognosticVariables{NF}(xvector,pvector,svector)
end