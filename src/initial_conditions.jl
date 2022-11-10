"""
Struct holding the so-called 'prognostic' variables.
'Prognostic' terminology is borrowed from cliamte science where it refers to 
any variables that are predicted via integration
"""
struct PrognosticVariables{NF<:AbstractFloat}
    xvector         ::AbstractVector{NF}       
    pvector         ::AbstractVector{NF}      
    svector         ::AbstractVector{NF}       
end


"""
    initialization = initial_conditions(M)
Setup the initial conditions for the MPD orbital dynamics"""
function initial_conditions(M::Model)


    @unpack NF = M.parameters
   
    # 1. Four- position
    @unpack r_initial,θ_initial, ϕ_initial = M.constants
    xvector = [0.0,r_initial,θ_initial, ϕ_initial] # Starting coordinates

    #1.1 Define some useful quantities
    @unpack a = M.parameters
    r,θ= xvector[2],xvector[3] 
    Δ = delta(r,a)
    Σ = sigma(r,θ,a)


    #To allow mutation of elements.
    #This is not needed in the actual covariant_metric() function
    #Why is this? Perhaps down to how SciMLSensitivity/DifferentialEquations.jl handles things?
    #Cleaner to just all g = covariant_metric(coords,a)
    xs = zeros(typeof(a),4,4)
    metric_covar = Zygote.Buffer(xs) # https://fluxml.ai/Zygote.jl/latest/utils/#Zygote.Buffer
    metric_covar[1,1] =  metric_g11(xvector,a) 
    metric_covar[2,2] =  metric_g22(xvector,a) 
    metric_covar[3,3] =  metric_g33(xvector,a) 
    metric_covar[4,4] =  metric_g44(xvector,a) 
    metric_covar[1,4] =  metric_g14(xvector,a) 
    metric_covar[4,1] =  metric_covar[1,4]
    g = copy(metric_covar)

   
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

    

    # #3. Four - spin
    @unpack Sθ,Sϕ = M.parameters
    @unpack s0 = M.constants 

    #Set the spatial components of the spin vector
    svector2 =  s0 * sin(Sθ) * cos(Sϕ)
    svector3 = -s0 *cos(Sθ)/r 
    svector4 =  s0*sin(Sθ)*sin(Sϕ)/(r*sin(θ)) 
    svector1 = -(svector2*pvector_covar[2] + svector3*pvector_covar[3]+svector4*pvector_covar[4])/pvector_covar[1] # Enforces the spin condition

    svector = [svector1,svector2,svector3,svector4]

    #conversion to NF happens here
    return PrognosticVariables{NF}(xvector,pvector,svector)

end