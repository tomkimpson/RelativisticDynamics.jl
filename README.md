# RelativisticDynamics

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tomkimpson.github.io/RelativisticDynamics.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tomkimpson.github.io/RelativisticDynamics.jl/dev/)
[![Build Status](https://github.com/tomkimpson/RelativisticDynamics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tomkimpson/RelativisticDynamics.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/tomkimpson/RelativisticDynamics.jl/branch/main/graph/badge.svg?token=LpfCcTaxFQ)](https://codecov.io/gh/tomkimpson/RelativisticDynamics.jl)

Welcome to RelativisticDynamics.jl!

The focus of this package is to numerically solve orbital equations for the Kerr spacetime.


Please see the [documentation](https://tomkimpson.github.io/RelativisticDynamics.jl/dev/)



## Spherical Photon Orbits

The spherical photon orbits are included as a test case or toy model to ensure the code machinery works as expected.

The general description of the equations of motion is described in [Teo 2003](https://ui.adsabs.harvard.edu/abs/2003GReGr..35.1909T/abstract) which have been reframed here via a Hamiltonian formalism c.f. [Pu et al. 2016](https://arxiv.org/abs/1601.02063)


## Relativistic spin orbital dynamics

The primary focus of this package is the calculation of the orbital motion of an astrophysical body, such as a pulsar, on a background Kerr spacetime. Going beyond the point particle approximation, we instead model an extended spinning body via the Mathisson-Papertrou-Dixon equations.

The code solves a set of ODEs numerically. These equations are based on the original works of Mathisson 1937, Papetrou 1951 and Dixon 1964. Consequently these equations are known as the MPD equations. More recent works can be found in Mashoon & Singh 2006, Singh 2005 and Singh, Wu & Sarty 2014.

Additional interesting discussion on the motion of extended bodies in GR can be found in Consta and Natatio, 2015



# Getting started 

```
solution,model = orbit();
```

