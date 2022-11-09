# Initial Conditions

We need to initialise 3 vectors $x^{\mu}, p^{\mu}, s^{\mu}$ for position, momentum and spin respectively. All this takes place in `src/universal_constants.jl` and `src/initial_conditions.jl`.


## Position $x^{\mu}$

The initial coordinates are specified straightforwardly. [Working in Boyer-Lindquist coordinates](https://en.wikipedia.org/wiki/Boyer%E2%80%93Lindquist_coordinates), we set $t=0$, $\theta = \pi/2$, $\phi=0.0$ and set the radial coordiante $r$ equal to the value of the semi-major axis specified in `src/system_parameters.jl`. 


## Momentum $p^{\mu}$
The specification of the initial momentum proceeds via two steps.

1. Map from the user-specified Keplerian orbital parameters to the conserved quantities
2. Define the momenta from the first-order equations of motion derived from the Hamilton-Jacobi function for the Kerr metric.


Both of these steps are well-described in [Schmidt 2002](https://arxiv.org/abs/gr-qc/0202090). 

To summarise the main steps:


For the Kerr spacetime we have 3 conserved quantities, the energy $E$, the angular momentum $L_z$ and the Carter constant $Q$.

Given the first order ODEs for $r$ and $\theta$ from the Kerr Hamiltonian, solve 

$$  \frac{dr}{d\lambda} = 0 ; \, \, \frac{d\theta}{d\lambda} = 0 $$

i.e. the turning points of the radial and polar motion. One can solve these equations to find $E,L,Q$ given $\alpha, e, \iota$. 


With the conserved quantities in hand, the 4-velocity is defined from the Kerr Hamiltonian,


$$\sigma \frac{dt}{d\lambda} = \frac{r^2 + a^2}{\Delta} P - a(aE\sin^2 \theta -L_z) $$





To summarise, 


