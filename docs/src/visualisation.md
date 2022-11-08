# Viewing your solution

With a solution in hand this can be quickly plotted using the usual [Plots.jl](https://docs.juliaplots.org/stable/) interface e.g.

```julia

using Plots
plot(solution,idxs=[1,2]) #to plot the time evolution of the first and second variables, in this case t and r 
```

Alternatively the user can use the inbuilt `PlotTrajectory` to plot the spatial path in usual`x,y,z` coords, in either 3D or 2D 

```julia


PlotTrajectory(solution,model,[1,2],"example_media/2D_example.png") #Plot a 2D example
PlotTrajectory(solution,model,[1,2,3],"example_media/3D_example.png") #Plot a 3D example


```