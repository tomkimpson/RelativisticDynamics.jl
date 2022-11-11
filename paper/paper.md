---
title: 'RelativisticDynamics.jl: Relativistic Spin-Orbital Dynamics in Julia'
tags:
  - Julia
  - astronomy
  - dynamics
  - gravity
  - general relativity
authors:
  - name: Tom Kimpson
    orcid: 0000-0002-6542-6032
    #equal-contrib: true
    affiliation: 1, 2 # (Multiple affiliations must be quoted)
  # - name: Tom Kimpson 2
  #   equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
  #   affiliation: 2
affiliations:
 - name: School of Physics, University of Melbourne, Parkville, VIC 3010, Australia 
   index: 1
 - name: Australian Research Council (ARC) Centre of Excellence for Gravitational Wave Discovery (OzGrav)
   index: 2
date: 14 October 2022
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---





<!-- A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience. -->

# Summary
Relativistic binaries composed of a millisecond pulsar (MSP) orbiting a much more massive ($\gtrsim 10^3 M_{\odot}$), spinning black hole (BH) are exceptional probes for investigating key questions of fundamental physics and astrophysics. Such systems are natural sources of gravitational waves (GWs) in the mHz regime, expected to be detectable by the next generation of space-based GW detectors such as LISA [@LISA]. The associated radio emission from the companion pulsar raises the possibility of an electromagnetic (EM) counterpart, enabling high precision multimessenger measurements to be made. The description of the orbital dynamics of these systems, and the influence on the resultant observed EM and GW signals, is non-trivial. A proper treatment of the spin-orbital dynamics can be derived from the conservation of the energy-momentum tensor
\begin{equation}\label{eq:conservation}
{T^{\mu \nu}}_{;\nu} = 0
\end{equation}
which when expanded into a set of infinite multipole moments leads to a description of the momentum vector $p^{\mu}$ and the spin tensor $s^{\mu \nu}$ 
\begin{equation}\label{eq:mpd1}
 \frac{Dp^{\mu}}{d \lambda} = -\frac{1}{2}{R^{\mu}}_{\nu \alpha \beta} u^{\nu} s^{\alpha \beta}
\end{equation}
\begin{equation}\label{eq:mpd2}
\frac{Ds^{\mu \nu}}{d \lambda} =p^{\mu}u^{\nu} - p^{\nu}u^{\mu}
\end{equation}
for affine parameter $\lambda$, 4-velocity $u^{\nu}$ and Riemann curvature tensor ${R^{\mu}}_{\nu \alpha \beta}$.  The system is closed by providing a spin supplementary condition
\begin{equation}\label{eq:mpd3}
s^{\mu \nu} p_{\nu} = 0
\end{equation}
Together, equations \ref{eq:mpd1} - \ref{eq:mpd3} form the Mathisson-Papetrou-Dixon (MPD) equations [@Mathisson1937;@Papapetrou1951; @Dixon1964], and describe the spin-orbital evolution in a fully consistent way that is applicable to strong field regimes. 


<!-- A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work -->
# Statement of need

`RelativisticDynamics.jl` is an open-source Julia package for relativistic spin-orbital dynamics in the gravitational strong field for a Kerr spacetime. Existing codes for modelling the dynamics of spinning objects like pulsars in the strong-field regime are generally lacking, since such systems occupy an intermediate regime that is generally overlooked. At the "low" end of this regime there are post-Newtonian or geodesic descriptions which neglect the influence of the pulsar spin on the underlying spacetime metric ("spin-curvature" coupling). At the "high" end there is the full Numerical Relativity (NR) solutions which are primarily applicable to two BHs with a mass ratio $\mathcal{O}(1)$, and are computationally intractable for these MSP systems which are observed over a large number of orbital cycles.  

`RelativisticDynamics.jl` aims to bridge this gap by providing a modern, fast code for accurate numerical evolution of spinning relativistic systems, via the MPD formalism. Julia is a modern language that solves the "two language problem", enabling fast dynamic typing and JIT compilation in conjunction with petaflop performance, comparable with numerical languages that are better known in the astrophysics community such as C or Fortran. As a modern language, it also provides a dedicated package manager and a large catalogue of _composable_ packages for scientific computing. This enables `RelativisticDynamics.jl` to easily leverage and interface with other scientific computing packages. The author and collaborators have used the general methods and mathematics described in this package for multiple research projects [e.g. @Kimpson2019; @Li2019; @Kimpson2020; @KimpsonAA] with a particular focus on the radio signals from spinning pulsar systems.  This package represents an attempt to create a documented, well-tested, open source resource for public use in this area, that can also be used as a computational playground for exploring techniques that could be applicable to more advanced numerical models.

In addition to providing a fast, modern package for strong field spin dynamics, `RelativisticDynamics.jl` has two additional important features from the perspective of modern relativistic astrophysics. Firstly, it is fully type flexible, being able to support arbitrary number formats. By making use of Julia's type-flexibility the model is written in such a way so as to be able to support hardware accelerated, low precision arithmetic and alternative rounding methods such as stochastic rounding. This enables rapid prototyping and exploration of reduced precision numerical techniques in astrophysics, an approach common in other numerical fields such as weather and climate modelling [e.g. @ECMWF]. Secondly, `RelativisticDynamcis.jl` is written to be fully differentiable via automatic differentiation. This enables the package to be used for differentiable physics applications in astrophysics, for example gravitational waveform modelling and parameter estimation or training neural networks based on the model. Automatic differentiation also provides a potential avenue for extension of the package to general (i.e. non-Kerr) spacetimes, whereby a user can specify the metric and the associated Christoffel symbols and Riemann tensors - which are simply linear combinations of the metric derivatives - are calculated automatically. 


Future potential extensions of this code include taking the dynamics beyond second order in the multipole expansion, and the inclusion of alternative spin conditions and spacetime metrics. The inclusion of a diagnostics tool for extracting gravitational waveforms in the time domain via a numerical kludge method would also be a worthwhile addition. Moreover, we have considered only bound dynamical systems - the ability to also explore hyberbolic systems would also be an interesting development.




![](../example_media/e01_stacked.pdf){width=50%}
![](../example_media/e08_stacked.pdf){width=50%}
\begin{figure}[!h]
\caption{Example orbital trajectories for a ms-pulsar with eccentricity $e=0.1$ (left panels), $e=0.8$ (right panels), orbiting a massive BH with extremal spin, $a=0.998$. The orbital motion is presented in the $x-y$ plane (top panels) and $x-z$ plane (bottom panels). The pulsar is initialised in the orbital plane with zero inclination. In the absence of spin-curvature coupling the particle would remain in the plane ($z=0$). Note the $z$-motion is on the scale of km, not gravitational radii. \label{fig:example}}
\end{figure}


# Acknowledgements

This work exploring the spin-evolution of relativistic systems via the the MPD equations was originally motivated through interesting discussions with Kinwah Wu. The port to a modern, precision-flexible model in Julia was heavily inspired by Milan Klöwer. Our thanks to both


# References





