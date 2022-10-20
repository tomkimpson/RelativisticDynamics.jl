# RelativisticDynamics.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tomkimpson.github.io/RelativisticDynamics.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tomkimpson.github.io/RelativisticDynamics.jl/dev/)
[![Build Status](https://github.com/tomkimpson/RelativisticDynamics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tomkimpson/RelativisticDynamics.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/tomkimpson/RelativisticDynamics.jl/branch/main/graph/badge.svg?token=LpfCcTaxFQ)](https://codecov.io/gh/tomkimpson/RelativisticDynamics.jl)

Welcome to RelativisticDynamics.jl!

The focus of this package is to numerically solve orbital equations for the Kerr spacetime.

Please see the [documentation](https://tomkimpson.github.io/RelativisticDynamics.jl/dev/).

## Relativistic spin orbital dynamics

The primary focus of this package is the calculation of the orbital motion of a spinning astrophysical body, such as a pulsar, on a background Kerr spacetime. 

The code solves a set of ODEs numerically. These equations are based on the original works of [Mathisson 1937](https://link.springer.com/article/10.1007/s10714-010-0939-y), [Papapetrou 1951](https://royalsocietypublishing.org/doi/10.1098/rspa.1951.0200) and [Dixon 1964](https://ui.adsabs.harvard.edu/abs/1964NCim...34..317D). Consequently these equations are known as the [MPD equations](https://en.wikipedia.org/wiki/Mathisson%E2%80%93Papapetrou%E2%80%93Dixon_equations). More recent works can be found in [Mashoon & Singh 2006](https://arxiv.org/abs/astro-ph/0608278), [Singh 2005](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.72.084033), [Singh, Wu & Sarty 2014](https://arxiv.org/abs/1403.7171) and [Li,Wu & Singh 2019](https://arxiv.org/abs/1902.03146). Additional interesting discussion on the motion of extended bodies in GR can be found in [Costa & Natário, 2015](https://arxiv.org/abs/1410.6443)






# Getting started 

```
using RelativisticDynamics
solution,model = orbit();
```
`solution` is the time evolution of the variables $x^{\mu}, p^{\mu}, s^{\mu}$. `model` is a struct containing a copy of all the settings that were used to generate that solution.

All user-defined parameters are set in `system_parameters.jl`. Default values can be overridden by passing to `orbit()` e.g.

```
solution,model = orbit(NF=Float32,model=:MPD)
```


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





# Useful literature

[Schmidt 2002](https://arxiv.org/abs/gr-qc/0202090)

[Barausse et al. 2007](https://arxiv.org/abs/0704.0138v2)