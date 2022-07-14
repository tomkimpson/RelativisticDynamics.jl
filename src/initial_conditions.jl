"""Struct holding the so-called 'prognostic' variables"""
struct PrognosticVariables{NF<:AbstractFloat}
    xvector         ::AbstractVector{NF}       
    pvector         ::AbstractVector{NF}      
    svector         ::AbstractVector{NF}       
end



"""Initialize a PrognosticVariables struct for an atmosphere at rest. No winds,
hence zero vorticity and divergence, but temperature, pressure and humidity are
initialised """
function initial_conditions(M::Model)

    @unpack NF,α,a,mPSR,mBH,Sθ,Sϕ = M.parameters
    @unpack E,L,Q,s0 = M.constants 



    # 4- position
    xvector = [0.0,α,pi/2.0,0.0] # By default the starting coordinates


    r,θ = xvector[2],xvector[3] #Get the initial r and theta 
    Δ = delta(r,a)
    Σ = sigma(r,θ,a)


    # Metric 
    metric_covar  = covariant_metric(r,θ,a)
    metric_contra = contravariant_metric(metric_covar,Δ*sin(θ)^2)

    (metric_contra)

    # 4 - momentum 
    T = (r^2 + a^2)*(E*(r^2 + a^2) -a*L)/Δ -a*(a*E*sin(θ)^2 - L)
    R = ((r^2 + a^2)*E -a*L)^2 - Δ*(r^2 + (L -a*E) + Q)
    Θ = Q - ((1.0 - E^2)*a^2 + L^2/sin(θ)^2)*cos(θ)^2                                                 # Note that this is a capital \Theta, not a lower case \theta
    Φ = a*(E*(r^2 + a^2) -a*L)/Δ -a*E + L/sin(θ)^2 


    tdot = T/Σ
    rdot = sqrt(R)/Σ
    θdot = sqrt(Θ)/Σ
    ϕdot = Φ/Σ
    m0 = mPSR/mBH # BH mass normalized to unity, so this is the pulsar mass in the new units

    pvector = m0*[tdot,rdot,θdot,ϕdot]

    pvector_covar = convert_to_covariant(metric_covar,pvector)


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