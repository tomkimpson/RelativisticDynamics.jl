using LaTeXStrings
using RecipesBase

"""
Plotting recipe for use with Plots.jl

"""
@recipe function f(r::CartesianTrajectory,dimensions=[1,2])

    #Default plotting settings

    size      -->  (600, 600)
    linecolor --> 1
    legend    --> :none


    #User can specify which dimensions to plot [1,2,3] --> [x,y,z]
    position_labels = [L"x (r_h)",L"y (r_h)",L"z (r_h)"]
    xlabel    --> position_labels[dimensions[1]]
    ylabel    --> position_labels[dimensions[2]]

    # return data
    data = [r.x,r.y,r.z]
    data[dimensions[1]], data[dimensions[2]]

    r.x, r.y

end 


# """
# Plotting recipe for use with Plots.jl
# """
# @recipe function f(r::CartesianResults)

#     seriestype --> :scatter
#     #markershape --> :auto        # if markershape is unset, make it :auto
#     markercolor :=  :blue  # force markercolor to be customcolor
#     #xrotation   --> 45           # if xrotation is unset, make it 45
#     #zrotation   --> 90           # if zrotation is unset, make it 90

#     # get the seriescolor passed by the user
#     #c = get(plotattributes, :seriescolor, :auto)

#     # return data
#     r.x, r.y

# end 






