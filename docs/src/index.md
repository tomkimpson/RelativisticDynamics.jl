```@meta
CurrentModule = RelativisticDynamics
```

# RelativisticDynamics.jl documentation

Welcome to the documentation for [RelativisticDynamics.jl](https://github.com/tomkimpson/RelativisticDynamics.jl), a relativistic orbital dynamics model for simulating spinning objects in curved spacetime.

## Overview
RelativisticDynamics.jl is a numerical model for determining the spin-orbital dynamics of a spinning body, such as a pulsar, on a background Kerr spacetime. The code solves a set of ODEs numerically. These equations are based on the original works of [Mathisson 1937](https://link.springer.com/article/10.1007/s10714-010-0939-y), [Papapetrou 1951](https://royalsocietypublishing.org/doi/10.1098/rspa.1951.0200) and [Dixon 1964](https://ui.adsabs.harvard.edu/abs/1964NCim...34..317D). Consequently these equations are known as the [MPD equations](https://en.wikipedia.org/wiki/Mathisson%E2%80%93Papapetrou%E2%80%93Dixon_equations). More recent works can be found in [Mashoon & Singh 2006](https://arxiv.org/abs/astro-ph/0608278), [Singh 2005](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.72.084033), [Singh, Wu & Sarty 2014](https://arxiv.org/abs/1403.7171) and [Li,Wu & Singh 2019](https://arxiv.org/abs/1902.03146). Additional interesting discussion on the motion of extended bodies in GR can be found in [Costa & Natário, 2015](https://arxiv.org/abs/1410.6443)



## Manual Outline

Please see the following pages of the documentation for more details    


- [How to run RelativisticDynamics.jl](how_to_run.md)
- [Initial conditions](IC.md)
- [MPD Equations](how_to_run.md)
- [Visualisation](visualisation.md)
- [Index of functions](how_to_run.md)


## Scope

The focus of SpeedyWeather.jl is to develop a global atmospheric model of intermediate complexity, that can run at various levels of precision (16, 32 and 64-bit) on different architectures (x86 and ARM, currently planned, GPUs probably in the future). Additionally, the model is written in an entirely number format-flexible way, such that any custom number format can be used and Julia will compile to the format automatically.


## Installation

SpeedyWeather.jl is registered in the Julia Registry. Open Julia's package manager from the REPL with `]`
and `add` the github repository to install SpeedyWeather.jl and all dependencies
```julia
(@v1.7) pkg> add SpeedyWeather
```
which will automatically install the latest release. However, you may want to install directly from the
main branch with
```julia
(@v1.7) pkg> add https://github.com/milankl/SpeedyWeather.jl#main
```
other branches than `#main` can be installed by adding `#branch_name` instead.




## Acknowledgements 
This work solving the MPD equations was originally motivated through interesting discussions with [Kinwah Wu](https://www.ucl.ac.uk/mssl/people/prof-kinwah-wu). The port to a modern, precision-flexible model in Julia was heavily inspired by [Milan Klöwer](https://github.com/milankl). Huge thanks to both.


Contributions are always welcome - just [open a pull request](https://github.com/tomkimpson/RelativisticDynamics.jl/pulls)


