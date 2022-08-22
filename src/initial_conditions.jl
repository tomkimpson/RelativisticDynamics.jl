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
    @unpack L,Q,s0,m0 = M.constants 


    # 1. Four- position
    xvector = [0.0,α,θ,ϕ] # By default the starting coordinates


    r,θ = xvector[2],xvector[3] #Get the initial r and theta 
    Δ = delta(r,a)
    Σ = sigma(r,θ,a)
    g  = covariant_metric(xvector,a)
    #metric_contra = contravariant_metric(xvector,a)


    # 2. Four - momentum 
    Pbar = (r^2+a^2) - a*L 
    Tbar = (r^2+a^2)*Pbar/Δ -a*(a*sin(θ)^2-L)
    Rbar = ((r^2+a^2)  -a*L)^2 -Δ*(r^2 + (L-a)^2+Q)
    θbar = Q - (L^2/sin(θ)^2)*cos(θ)^2
    ϕbar = a*Pbar/Δ -a + L/sin(θ)^2

    #These are 4 velocities from Schmidt 2002
    tdot = -Tbar/Σ
    rdot = -sqrt(Rbar)/Σ
    θdot = -sqrt(θbar)/Σ
    ϕdot = -ϕbar/Σ

    pvector = m0*[tdot,rdot,θdot,ϕdot]
    pvector_covar = convert_to_covariant(g,pvector)




    #3. Four - spin
    svector = [0.0,0.0,0.0,0.0]

    #Set the spatial components of the spin vector
    svector[2] = s0 * sin(Sθ) * cos(Sϕ)
    svector[3] = -s0 *cos(Sθ)/r 
    svector[4] = s0*sin(Sθ)*sin(Sϕ)/(r*sin(θ)) 
    svector[1] = -(svector[2]*pvector_covar[2] + svector[3]*pvector_covar[3]+svector[4]*pvector_covar[4])/pvector_covar[1] # Enforces the spin condition


    # conversion to NF happens here
    return PrognosticVariables{NF}(xvector,pvector,svector)
end