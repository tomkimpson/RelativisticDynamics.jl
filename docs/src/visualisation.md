# Visualisation of the solution
It is often desireable, as a sanity check, to plot the orbital solution. Once the `orbit()` call completes, the output can be quickly plotted using the usual [Plots.jl](https://docs.juliaplots.org/stable/) interface e.g.

```julia
using Plots
using RelativisticDynamics

solution,model = orbit()
plot(solution,model.constants.a) 
```
By default the plotting recipe plots the $x-y$ orbital trajectory and up-samples the numerical solution by a factor of $10$. The user can define their own up-sampling factor, and pass any two of the integration variables as arguments, as well as the Cartesian coordinates ($x,y,z$). For example:

```julia
plot(solution,model.constants.a,upsample=2, vars_to_plot = [:t,:r])     # Plot a timeseries of the r-coordinate, upsampled by a factor of 2
plot(solution,model.constants.a,upsample=100, vars_to_plot = [:sθ,:sϕ]) # Plot the θ-ϕ components of the spin vector against each other, upsampled by a factor of 100
plot(solution,model.constants.a,vars_to_plot = [:x,:z]) # Plot the x-z orbital trajectory

```

