# How to run RelativisticDynamics.jl

The simplest way to run RelativisticDynamics.jl with default parameters is

```julia
using RelativisticDynamics
solution,model = orbit()
```

The `orbit()` function returns two objects. The first, `solution` holds the evolution of the position, momentum and spin vectors. The second, `model`, holds a copy of all the parameters and settings used to generate the solution, e.g. what was the BH spin?

All default parameters can be found in `src/system_parameters.jl`. Passing a keyword argument to `orbit()` overrides the defaults e.g.

```jldoctest
using RelativisticDynamics
solution,model = orbit(e=0.6,a=0.99);
model.constants.a
# output
0.99
```
generates the solution for a system with an eccentricity = 0.6, around a BH with an extremal spin. 

If provided, the number format has to be the first argument, all other arguments are keyword arguments. e.g. 

```jldoctest
using RelativisticDynamics
solution,model = orbit(Float32,e=0.6,a=0.99);
model.parameters.NF
# output
Float32
```


Please see `notebooks/demo.ipynb` for some worked examples using `RelativisticDynamics.jl`, including the application of autodiff methods.


