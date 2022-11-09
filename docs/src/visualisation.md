# Visualisation of the solution
It is often desireable, as a sanity check, to plot the orbital solution. Once the `orbit()` call completes, the output can be quickly plotted using the usual [Plots.jl](https://docs.juliaplots.org/stable/) interface e.g.

```julia
using Plots
using RelativisticDynamics

solution,model = orbit()  
plot(solution,idxs=[1,2]) # to plot the time evolution of the first and second variables, in this case t and r 
```

Alternatively the user can use the inbuilt `PlotTrajectory()` function to plot the spatial path in usual Cartesian `x,y,z` coords, in either 3D or 2D. Here the indexes 1,2,3 refer to the `x`,`y`,`z` coordinates respectively. 

```julia
PlotTrajectory(solution,model,[1,2],"example_media/2D_example.png")   # Plot a 2D example, save to example_media/2D_example.png
PlotTrajectory(solution,model,[1,2,3],"example_media/3D_example.png") # Plot a 3D example, save to example_media/2D_example.png
```


One can also create a stacked plot of the trajectory in the $x-y$ and $x-z$ planes:
```julia
StackedPlot(solution,model,"../example_media/e08_stacked.pdf")
```
