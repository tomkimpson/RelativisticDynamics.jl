# Initial Conditions

We need to initialise 3 vectors $x^{\mu}, p^{\mu}, s^{\mu}$ for position, momentum and spin respectively. All this takes place in `src/universal_constants.jl` and `src/initial_conditions.jl`.


## Position, $x^{\mu}$

The initial coordinates are specified straightforwardly. [Working in Boyer-Lindquist coordinates](https://en.wikipedia.org/wiki/Boyer%E2%80%93Lindquist_coordinates), we set $t=0$, $\theta = \pi/2$, $\phi=0.0$ and set the radial coordiante $r$ equal to the value of the semi-major axis $\alpha$ specified in `src/system_parameters.jl`. 


## Momentum, $p^{\mu}$
The specification of the initial momentum proceeds via two steps.

1. Map from the user-specified Keplerian orbital parameters to the conserved quantities
2. Define the momenta from the first-order equations of motion derived from the Hamilton-Jacobi function for the Kerr metric.


Both of these steps are well-described in [Schmidt 2002](https://arxiv.org/abs/gr-qc/0202090). The general overview is as follows:


For the Kerr spacetime we have 3 conserved quantities, the energy $E$, the angular momentum $L_z$ and the Carter constant $Q$. We want to be able to determine the value of these quantities, given the Keplerian orbital parameters, $\alpha, e, \iota$. 

Given the first order ODEs for $r$ and $\theta$ from the Kerr Hamiltonian, solve 

$$\frac{dr}{d\lambda} = 0 ; \, \, \frac{d\theta}{d\lambda} = 0$$
i.e. the turning points of the radial and polar motion. One can solve these equations to find $E,L,Q$ given $\alpha, e, \iota$. 

With the conserved quantities in hand, the 4-velocity is defined from the Kerr Hamiltonian,

$$\sigma \frac{dt}{d\lambda} = \frac{r^2 + a^2}{\Delta} P - a(aE\sin^2 \theta -L_z)$$

$$\sigma \frac{dr}{d\lambda} = \pm \sqrt{R}$$

$$\sigma \frac{d\theta}{d\lambda} = \pm \sqrt{\Theta}$$

$$\sigma \frac{d\phi}{d\lambda} = \frac{a}{\Delta} - aE + \frac{L_z}{\sin^2 \theta}$$
where $R,\Theta$ and $P$ are again given in [Schmidt 2002](https://arxiv.org/abs/gr-qc/0202090). We always take the positive square root for the initialisation, such that the initial motion of the pulsar is "outwards and upwards" (increasing $r$ and $\theta$). 

The covariant 4-velocity can then be translated into a contravariant 4-momentum as
$$p^{\alpha} = m g^{\alpha \beta} u_{\beta}$$
for metric $g^{\alpha \beta}$.
## Spin, $s^{\mu}$

In order to determine the spin vector we must first specify the moment of inertia of the pulsar. We model the pulsar as a solid sphere such that

$$I = \frac{2}{5} m_{\rm PSR} r_{\rm PSR}^2$$
The angular momentum/spin magnitude is, 

$$s_0 = 2 \pi I / P_{\rm PSR}$$
where $P_{\rm PSR}$ is the spin period of the pulsar. The spatial components of the spin vector are then,

$$s^r = s_0 \sin(S_{\theta}) \cos(S_{\phi})$$

$$s^{\theta} = -s_0 \cos(S_{\theta})/r$$

$$s^{\phi} = s_0 \sin(S_{\theta}) \sin(S_{\phi})/r \sin(\theta)$$
where $S_{\theta, \phi}$ are the latitude and azimuthal angles of the spin axis, see e.g. [Mashhoon & Singh, 2006](https://arxiv.org/abs/astro-ph/0608278). The temporal component $s^{t}$ is enforced by the spin condition. Throughout this package we take the Tulczyjew-Dixon (TD) condition (see e.g. [Costa & Natário, 2015](https://arxiv.org/abs/1410.6443) for discussion of TD condition and other options)
 $$s^{\mu}p_{\mu}  = 0 $$











