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
function BoyerLindquistPlot(solution,M)

    @unpack a = M.parameters

    r = solution[2,:]
    θ = solution[3,:]
    ϕ = solution[4,:]


    w = sqrt.(r.^2 .+ a^2) 
    x = w .* sin.(θ) .* cos.(ϕ)
    y = w .* sin.(θ) .* sin.(ϕ)
    z = r .* cos.(θ)




    plt = plot3d(1,title = "TEST TITLE", marker=2,markercolor="black")

     n = length(x)
     anim = @animate for i in 1:n
        push!(plt,x[i],y[i],z[i])
     end 
     gif(anim, "yt2.gif", fps = 30)





    #Set the backend 
    #gr()
    
    #plotly(lw=3)
    #plot(x,y)
    #scatter(x,y,z)
    #plot(solution,vars=[1])

end 





@userplot CirclePlot
@recipe function f(cp::CirclePlot)
    x, y, i = cp.args
    n = length(x)
    inds = circshift(1:n, 1 - i)
    linewidth --> range(0, 10, length = n)
    seriesalpha --> range(0, 1, length = n)
    aspect_ratio --> 1
    label --> false
    x[inds], y[inds]
end



function PlotSphericalPhotonOrbit(solution,model)

    @unpack a = model.parameters
    @unpack u0 = model.constants


    T = range(first(solution.t),last(solution.t),length=1000)

    p = solution(T)
    r = p[2,:]
    χ = p[3,:] # reparameterised
    ϕ = p[4,:]

    #Convert back to θ
    θ = acos.(sqrt(u0)*cos.(χ)) 

    w = sqrt.(r.^2 .+ a^2) 
    x = w .* sin.(θ) .* cos.(ϕ)
    y = w .* sin.(θ) .* sin.(ϕ)
    z = r .* cos.(θ)


    #n = 150
    #t = range(0, 2π, length = n)
    #x = sin.(t)
    #y = cos.(t)
    
    # anim = @animate for i ∈ 1:length(x)
    #     circleplot(x, y, i,line_z = 1:length(x))
    # end
    # gif(anim, "anim_fps15.gif", fps = 15)


    plt = plot3d(1,title = "TEST TITLE", marker=2,markercolor="black")


     n = length(x)
     anim = @animate for i in 1:n
        push!(plt,x[i],y[i],z[i])
     end 
     gif(anim, "yt.gif", fps = 30)



    # anim = @animate for i ∈ 1:length(x)
    #     circleplot(x, y, i)
    # end

    #plot(x,y,z)
#     df=1
#     anim = @animate for i = 1:df:length(x)
#        plot(x[1:i], y[1:i],z[1:i], legend=false)
#    end
     
#     gif(anim, "tutorial_anim_fps30.gif", fps = 30)


    #plotlyjs()

    #scatter(x,y,z)

end 