


using Plots 
using LaTeXStrings
using Printf

"""
Plot the 3D trajectory of a body. Assumes coordinates are Boyer Lindquist. Saves a low resolution
PNG figure to disk
"""
function PlotTrajectory(solution,model,saveit)

    @unpack a,model = model.parameters    #Get the BH spin parameter 
    


    #Interpolate to higher resolution for smooth plotting   
    interpolation_factor = 10 
    T = range(first(solution.t),last(solution.t),length=length(solution.t)*interpolation_factor)
    p = solution(T)

    # Extract relevant data from the interpolated solution 
    r = p[2,:]
    θ = p[3,:] 
    ϕ = p[4,:]

    # Boyer lindquist to Cartesian 
    w = sqrt.(r.^2 .+ a^2) 
    x = w .* sin.(θ) .* cos.(ϕ)
    y = w .* sin.(θ) .* sin.(ϕ)
    z = r .* cos.(θ)


    # Define the horizon surface 
    n = 100
    u = range(-π, π; length = n)
    v = range(0, π; length = n)
    rH = 1.0 + sqrt(1.0 - a^2)
    xH = rH*cos.(u) * sin.(v)'
    yH = rH*sin.(u) * sin.(v)'
    zH = rH*ones(n) * cos.(v)'

    #Define plot limits
    #m = max()
    m = max(maximum(abs.(x)),maximum(abs.(y)),maximum(abs.(z)))
    #m=10

    # Plot it 
    

    #gr()
    #pyplot()  # Set the backend

    #plotlyjs()
    println(model)
    title = "Spherical photon orbits with a = $(@sprintf("%.2f", a))"
    plot(x,y,z,
              xaxis=(L"x (r_h)",(-m,m)),yaxis=(L"y (r_h)",(-m,m)),zaxis=(L"z (r_h)",(-m,m)),
              legend=false,
              title = title,
              size = (1200, 800))

 

    # # # Singularity 
    # pobject = plot!([0.0],[0.0],[0.0],marker=10)



    # pobject = wireframe!(xH, yH, zH) #or surface! 

    # if saveit
    #     fout = "example_media/spherical_photon_orbits_a_$a" * ".png"
    #     savefig(fout)
    # end


end 


# """
# Create a 3D gif of the trajectory of a body. Assumes coordinates are Boyer Lindquist.
# """
# function AnimateTrajectory(solution,model)

#     @unpack a = model.parameters
#     @unpack rH = model.constants

#     println("This is animate trajectory")

#     #Interpolate to higher resolution for smooth plotting   
#     interpolation_factor = 10 
#     T = range(first(solution.t),last(solution.t),length=length(solution.t)*interpolation_factor)
#     p = solution(T)

#     # Extract relevant data from the interpolated solution 
#     r = p[2,:]
#     θ = p[3,:] 
#     ϕ = p[4,:]

#     # Boyer lindquist to Cartesian 
#     w = sqrt.(r.^2 .+ a^2) 
#     x = w .* sin.(θ) .* cos.(ϕ)
#     y = w .* sin.(θ) .* sin.(ϕ)
#     z = r .* cos.(θ)



#     # Plot it 
#     title = "Spherical photon orbits with a = $(@sprintf("%.2f", a))"

#     plt = plot3d(1,
#                 xaxis=(L"x (r_h)",(-3,3)),yaxis=(L"y (r_h)",(-3,3)),zaxis=(L"z (r_h)",(-3,3)),
#                 legend=false,
#                 title = title,
#                 size = (1200, 800))


#      n = length(x)
#      anim = @animate for i in 1:n
#         push!(plt,x[i],y[i],z[i])
#      end 
#     fout = "example_media/spherical_photon_orbits_a_$a" * ".gif"
#      gif(anim, fout, fps = 30)


# end 







