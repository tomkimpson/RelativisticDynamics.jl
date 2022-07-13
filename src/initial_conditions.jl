"""Struct holding the so-called 'prognostic' variables"""
struct PrognosticVariables{NF<:AbstractFloat}
    xvector         ::AbstractVector{NF}       
    pvector         ::AbstractVector{NF}       

   # metric_covar    :: Array{NF,2}
   #metric_contra   :: Array{NF,2}

end



"""Initialize a PrognosticVariables struct for an atmosphere at rest. No winds,
hence zero vorticity and divergence, but temperature, pressure and humidity are
initialised """
function initial_conditions(M::Model)

    @unpack NF,α,a,mPSR,mBH = M.parameters
    @unpack E,L,Q = M.constants 



    # 4- position
    xvector = [0.0,α,pi/2.0,0.0] # By default the starting coordinates


    r,θ = xvector[2],xvector[3] #Get the initial r and theta 
    Δ = delta(r,a)
    Σ = sigma(r,θ,a)


    # Metric 
    metric_covar  = covariant_metric(r,θ,a)
    metric_contra = contravariant_metric(r,θ,a)


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



    # 4 - spin
    #update_metric!(M,r,θ)
    # metric = initialize_metric(M)
    # println(metric.metric_covar)
    # println(metric.metric_contra)



    # conversion to NF happens here
    return PrognosticVariables{NF}(xvector,pvector)
end