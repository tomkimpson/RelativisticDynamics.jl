# How to run RelativisticDynamics.jl

The simplest way to run RelativisticDynamics.jl with default parameters is

```julia
using RelativisticDynamics
solution,model = orbit()
```

The `orbit()` funciton returns two objects. The first, `solution` holds the evolution of the position, momentum and spin vectors. The second, `model`, holds a copy of all the parameters and settings used to generate the solution (e.g. what was the BH spin?).

All default parameters can be found in `src/default_parameters.jl`. Passing a keyword argument to `orbit()` overrides the defaults e.g.

```julia
run_speedy(e=0.6,a=0.99)
```
would generate the solution for a system with an eccentricity = 0.6, around a BH with an extramal spin. 

If provided, the number format has to be the first argument, all other arguments are keyword arguments. e.g. 
```julia
run_speedy(Float32,e=0.6,a=0.99)
```

Please see `notebooks/demo.ipynb` for some worked examples using `RelativisticDynamics.jl`, including the application of autodiff methods.



