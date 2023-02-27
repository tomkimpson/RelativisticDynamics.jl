# Mathisson-Papetrou-Dixon Equations

The MPD equations are derived from the conservation of the energy-momentum tensor:

$${T^{\mu \nu}}_{;\nu} = 0$$
Taking the multipole expansion leads to a description of the momentum vector $p^{\mu}$ (0th moment) and the spin tensor $s^{\mu \nu}$ (dipole moment)

$$\frac{Dp^{\mu}}{d \lambda} = -\frac{1}{2}{R^{\mu}}_{\nu \alpha \beta} u^{\nu} s^{\alpha \beta}$$
$$\frac{Ds^{\mu \nu}}{d \lambda} =p^{\mu}u^{\nu} - p^{\nu}u^{\mu}$$
for affine parameter $\lambda$, 4-velocity $u^{\nu}$ and Riemann curvature tensor ${R^{\mu}}_{\nu \alpha \beta}$, with $\frac{D}{d\lambda}$ denoting a covariant derivative. We ignore the higher order moments since for pulsar mass ($m$) << BH mass ($M$) and pulsar radius << BH radius, the motion is dominated by the lowest order moments. 

The system is closed by providing a spin supplementary condition (SSC). Choosing an SSC is equivalent to choosing the observer who measures the centre of mass of the spinning body - in GR the centre of mass is observer dependent! Multiple SSCs exist, and different choices of SSC will naturally lead to different solutions for the system. However, for the Kerr spacetime and within the pole-dipole approximation (see below) all these solutions for the centre of mass (the "minimal worldtube") are contained within the convex hull of the body's worldtube (see [Costa & Natário, 2015](https://arxiv.org/abs/1410.6443) for an thorough discussion of the spin conditions). Throughout this package we take the Tulczyjew-Dixon (TD) condition

$$s^{\mu \nu} p_{\nu} = 0$$
which is equivalent to choosing the centre of mass to be the one measured in the zero-3 momentum frame, which has the advantage of specifying a unique worldline 


In the extreme mass ratio limit m << M, the pulsar Möller radius is much less than the gravitational lengthscale. This means that the pole-dipole interaction is much stronger than the dipole-dipole interaction. Within this approximation, the MPD equations reduced to a set of ODEs (see e.g. [Mashoon & Singh, 2006](https://arxiv.org/abs/astro-ph/0608278), [Singh, Wu & Sarty, 2014](https://arxiv.org/abs/1403.7171))

$$\frac{dp^{\alpha}}{d\lambda} = - \Gamma_{\mu\nu}^{\alpha} p^{\mu}u^{\nu} +  \left( \frac{1}{2m} {R^{\alpha}}_{\beta \rho \sigma} \epsilon^{\rho \sigma}_{\quad \mu \nu} s^{\mu} p^{\nu} u^{\beta}\right) \ ,$$

$$\frac{ds^{\alpha}}{d \lambda} = - \Gamma^{\alpha}_{\mu \nu} s^{\mu}u^{\nu} + \left(\frac{1}{2m^3}R_{\gamma \beta \rho \sigma} \epsilon^{\rho \sigma}_{\quad \mu \nu} s^{\mu} p^{\nu} s^{\gamma} u^{\beta}\right)p^{\alpha} \ ,$$

$$\frac{dx^{\alpha}}{d\lambda} = -\frac{p^{\delta}u_{\delta}}{m^2} \left[ p^{\alpha} + \frac{1}{2} \frac{ (s^{\alpha \beta} R_{\beta \gamma \mu \nu} p^{\gamma} s^{\mu \nu})}{m^2 + (R_{\mu \nu \rho \sigma} s^{\mu \nu} s^{\rho \sigma}/4)}\right]$$
Whilst $\lambda$ has the freedom to be any affine parameter, we take it to be the proper time $\tau$ such that $g_{\mu \nu}u^{\mu} u^{\nu}=-1$. These equations are generally covariant - we take the astrophysically motivated [Kerr metric](https://en.wikipedia.org/wiki/Kerr_metric) for a spinning BH and work in Boyer-Lindquist coordinates.
