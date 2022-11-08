# RelativisticDynamics.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tomkimpson.github.io/RelativisticDynamics.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tomkimpson.github.io/RelativisticDynamics.jl/dev/)
[![Build Status](https://github.com/tomkimpson/RelativisticDynamics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tomkimpson/RelativisticDynamics.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/tomkimpson/RelativisticDynamics.jl/branch/main/graph/badge.svg?token=LpfCcTaxFQ)](https://codecov.io/gh/tomkimpson/RelativisticDynamics.jl)

Welcome to RelativisticDynamics.jl!

The focus of this package is to numerically solve orbital equations for the Kerr spacetime.

Please see the [documentation](https://tomkimpson.github.io/RelativisticDynamics.jl/dev/).

## Relativistic spin orbital dynamics






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





