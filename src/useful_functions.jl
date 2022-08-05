"""
Δ = delta(r,a)
The well-known delta function of the Kerr metric
"""
function delta(r,a)
return r^2 -2.0*r + a^2
end 

"""
Σ = sigma(r,θ,a)
The well-known sigma function of the Kerr metric
"""
function sigma(r,θ,a)
return r^2 + a^2 * cos(θ)^2
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
return r*(r-2.0) + (zminus^2)/(1.0 - zminus^2) *Δ
end


"""
d = mapping_d(r,a,zminus)
Mapping function `d` used when converting from Keplerian orbital parameters to constants of motion
"""
function mapping_d(r,a,zminus)
Δ = delta(r,a)
return (r^2 +a^2 * zminus^2)*Δ
end



"""
Convert a vector from contravariant form to convariant form using the metric 
"""
function convert_to_covariant(metric,vector)


    vector_covar = vector
    vector_covar[1] = metric[1,1]*vector[1] + metric[1,2]*vector[2] + metric[1,3]*vector[3] + metric[1,4]*vector[4]
    vector_covar[2] = metric[2,1]*vector[1] + metric[2,2]*vector[2] + metric[2,3]*vector[3] + metric[2,4]*vector[4]
    vector_covar[3] = metric[3,1]*vector[1] + metric[3,2]*vector[2] + metric[3,3]*vector[3] + metric[3,4]*vector[4]
    vector_covar[4] = metric[4,1]*vector[1] + metric[4,2]*vector[2] + metric[4,3]*vector[3] + metric[4,4]*vector[4]

    return vector_covar  

end 


using Plots 
using LaTeXStrings
using Printf

"""
Plot the 3D trajectory of a body. Assumes coordinates are Boyer Lindquist. Saves a low resolution
PNG figure to disk
"""
function PlotTrajectory(solution,model)

    @unpack a = model.parameters
    @unpack rH = model.constants


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



    # Plot it 
    pyplot()  # Set the backend

    title = "Spherical photon orbits with a = $(@sprintf("%.2f", a))"
    plot(x,y,z,
         xaxis=(L"x (r_h)",(-3,3)),yaxis=(L"y (r_h)",(-3,3)),zaxis=(L"z (r_h)",(-3,3)),
         legend=false,
         title = title,
         size = (1200, 800))

    # Singularity 
    plot!([0.0],[0.0],[0.0],marker=10)


    # Horizon surface 
    
    n = 100
    u = range(-π, π; length = n)
    v = range(0, π; length = n)
    xH = rH*cos.(u) * sin.(v)'
    yH = rH*sin.(u) * sin.(v)'
    zH = rH*ones(n) * cos.(v)'

    wireframe!(xH, yH, zH) #or surface! 

    fout = "example_media/spherical_photon_orbits_a_$a" * ".png"
    savefig(fout)
end 


"""
Create a 3D gif of the trajectory of a body. Assumes coordinates are Boyer Lindquist.
"""
function AnimateTrajectory(solution,model)

    @unpack a = model.parameters
    @unpack rH = model.constants

    println("This is animate trajectory")

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



    # Plot it 
    title = "Spherical photon orbits with a = $(@sprintf("%.2f", a))"

    plt = plot3d(1,
                xaxis=(L"x (r_h)",(-3,3)),yaxis=(L"y (r_h)",(-3,3)),zaxis=(L"z (r_h)",(-3,3)),
                legend=false,
                title = title,
                size = (1200, 800))


     n = length(x)
     anim = @animate for i in 1:n
        push!(plt,x[i],y[i],z[i])
     end 
    fout = "example_media/spherical_photon_orbits_a_$a" * ".gif"
     gif(anim, fout, fps = 30)


end 





















# function BoyerLindquistPlot(solution,M)

#     @unpack a = M.parameters

#     r = solution[2,:]
#     θ = solution[3,:]
#     ϕ = solution[4,:]


#     w = sqrt.(r.^2 .+ a^2) 
#     x = w .* sin.(θ) .* cos.(ϕ)
#     y = w .* sin.(θ) .* sin.(ϕ)
#     z = r .* cos.(θ)



#     plot(x,y,z)

#     # plt = plot3d(1,xaxis=("x",(-3,3)),yaxis=("y",(-3,3)),zaxis=("z",(-3,3)),title = "TEST TITLE", marker=2,markercolor="black")

#     #  n = length(x)
#     #  anim = @animate for i in 1:n
#     #     push!(plt,x[i],y[i],z[i])
#     #  end 
#     #  gif(anim, "example_media/spherical_photon_orbit.gif", fps = 30)





#     #Set the backend 
#     #gr()
    
#     #plotly(lw=3)
#     #plot(x,y)
#     #scatter(x,y,z)
#     #plot(solution,vars=[1])

# end 





# function PlotSphericalPhotonOrbit(solution,model)

#     @unpack a = model.parameters
#     @unpack u0 = model.constants


#     nsteps = 1000
#     T = range(first(solution.t),last(solution.t),length=1000)

#     println("Plotting over the time range ", first(solution.t)," ",last(solution.t), " with ", nsteps, " steps")

#     p = solution(T)
#     r = p[2,:]
#     θ = p[3,:] 
#     ϕ = p[4,:]

#     w = sqrt.(r.^2 .+ a^2) 
#     x = w .* sin.(θ) .* cos.(ϕ)
#     y = w .* sin.(θ) .* sin.(ϕ)
#     z = r .* cos.(θ)


#     #plotlyjs()
#     pyplot()
#     plot(x,y,z)


#     #Animation

#     # plt = plot3d(1,
#     #             xaxis=("x",(-3,3)),yaxis=("y",(-3,3)),zaxis=("z",(-3,3)),
#     #             title = "TEST TITLE", marker=2,markercolor="black")


#     #  n = length(x)
#     #  anim = @animate for i in 1:n
#     #     push!(plt,x[i],y[i],z[i])
#     #  end 
#     #  gif(anim, "yt.gif", fps = 30)



#     # anim = @animate for i ∈ 1:length(x)
#     #     circleplot(x, y, i)
#     # end

#     #plot(x,y,z)
# #     df=1
# #     anim = @animate for i = 1:df:length(x)
# #        plot(x[1:i], y[1:i],z[1:i], legend=false)
# #    end
     
# #     gif(anim, "tutorial_anim_fps30.gif", fps = 30)


#     #plotlyjs()

#     #scatter(x,y,z)

# end 