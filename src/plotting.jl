


using Plots 
using LaTeXStrings
using Printf

"""
Plot trajectory of a body. Assumes coordinates are Boyer Lindquist. 
Plots in either 2D or 3D depending on specification of dimensions.
Saves a low resolution PNG figure to disk in example_media/
"""
function PlotTrajectory(solution,model,dimensions=[1,2,3],savepath="")

    @unpack a = model.parameters    #Get the BH spin parameter 
    
    println("Plotting the solution generated with the following user defined parameters")
    display(model.parameters)
    println("-------------------------------")

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
    position = [x,y,z]
    position_labels = [L"x (r_h)",L"y (r_h)",L"z (r_h)"]


    #Setup plotting env
    if length(dimensions) == 3
        
        plot(x,y,z,
            legend=false,
            xlabel=position_labels[1],
            ylabel=position_labels[2],
            zlabel=position_labels[3],
            camera = (25, 30),
            size = (1000, 600))


        xBH = 0:0; yBH = 0:0; zBH = 0:0
        scatter!(xBH, yBH,zBH,markercolor="red",markersize=5) 

    elseif length(dimensions) == 2
        idx1,idx2 = dimensions 
        plt = plot(position[idx1],position[idx2],
             xlabel=position_labels[idx1],
             ylabel=position_labels[idx2],
             legend=false,
             size = (600, 600)
             )

        xBH = 0:0; yBH = 0:0
        plt = scatter!(xBH, yBH,markercolor="red",markersize=5)

    else
        println("Those dimensions are not defined")
        return

    end

    display(plt)
    if ~isempty(savepath)
        println("Saving figure to: ", savepath)
        savefig(savepath)
    end



end 


"""
Plot the x-y and x-z trajectory of a body on two separate subplots
"""
function StackedPlot(solution,model,savepath="")


    #Interpolate to higher resolution for smooth plotting   
    @unpack a = model.parameters    #Get the BH spin parameter 

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
    

    #Convert to km in the z-direction
    mBH = model.parameters.mBH
    Msolar = 2e30
    c = 3e8
    G = 6.67e-11
    factor = G*mBH*Msolar/c^2
    z_km = z * factor/1e3


    #x-y 
    plot(x,y,layout=grid(2,1, heights=(0.8,0.2)), size=(450,600),legend=false,link = :x)
    plot!(ylabel=L"y (r_h)",subplot=1)



    xBH = 0:0; yBH = 0:0
    plt = scatter!(xBH, yBH,markercolor="red",markersize=5,subplot=1)

    #x-z
    plot!(x,z_km,subplot=2,legend=false)
    plot!(xlabel=L"x \rm (r_h)",subplot=2)
    plot!(ylabel=L"z (km)",subplot=2)


    display(plt)
    if ~isempty(savepath)
        println("Saving figure to: ", savepath)
        savefig(savepath)
    end



end 

