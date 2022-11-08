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
    orcid: 0000-0000-0000-0000
    equal-contrib: true
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
Relativistic binaries composed of a millisecond pulsar (MSP) orbiting a much more massive BH ($\gtrsim 10^3 M_{\odot}$) are exceptional probes for investigating key questions of fundamental physics and astrophysics. Such systems are natural sources of gravitational waves (GWs)in the mHz regime, expected to be detectable by the next generation of space-based GW detectors such as LISA \cite{LISA}. The associated radio emission from the companion pulsar raises the possibility of an electromagnetic counterpart, enabling high precision multimessenger measurements to be made. The description of the orbital dynamics of these systems, and the influence on the resultant EM and GW observed signal, is non-trivial. Whilst models commonly treat the motion via a post-Newtonian or geodesic description, such an approach neglects the non-linear influence of the pulsar spin on the underlying spacetime metric. A proper treatment of the spin-orbital dynamics can be derived from the conservation of the energy-momentum tensor
\begin{equation}\label{eq:conservation}
{T^{\mu \nu}}_{;\nu} = 0
\end{equation}
which when expanded into a set of infinite multipole moments leads to a description of the momentum vector $p^{\mu}$ and the spin tensor $s^{\mu \nu}$ 
\begin{equation}\label{eq:mpd1}
 \frac{Dp^{\mu}}{d \lambda} = -\frac{1}{2}R^{\mu}_{\nu \alpha \beta} u^{\nu} s^{\alpha \beta}
\end{equation}
\begin{equation}\label{eq:mpd2}
\frac{Ds^{\mu \nu}}{d \lambda} =p^{\mu}u^{\nu} - p^{\nu}u^{\mu}
\end{equation}
The system is closed by providing a spin supplementary condition
\begin{equation}\label{eq:mp3}
s^{\mu \nu} p_{\nu} = 0
\end{equation}
Together, equations \autoref{eq:md1} - \autoref{eq:md3} form the Mathisson-Papetrou-Dixon (MPD) equations, and describe the spin-orbital evolution in a fully consistent way that is applicable to strong field regimes. 

\ref{eq:mpd1}



<!-- A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work -->
# Statement of need

`RelativisticDynamics.jl` is an open-source Julia package for relativistic spin-orbital dynamics in the gravitational strong field. Existing codes for modelling the dynamics of spinning objects like pulsars in the strong-field regime are generally lacking, since such systems occupy an intermediate regime that is generally overlooked. At the "low" end of the gravitational there are post-Newtonian or geodesic descriptions which neglect non-linear spin contributions, whilst at the ""high" end there is the full are full Numerical Relativity (NR) solutions which are primarily applicable to two BHs with a mass ratio $\mathcal{O}(1)$, and are computationally intractable for these MSP systems which are observed over a large number of orbital cycles. 

`RelativisticDynamcis.jl` aims to bridge this gap by providing a modern, fast code for accurate numerical evolution of spinning relativistic systems, via the MPD formalism. Julia is a modern language that solves the "two language problem", enabling fast dynamic typing and JIT compilation in conjunction with petaflop performance, comparable with better known numerical languages such as C or Fortran. As a modern language, it also provides a dedicated package manager and a large catalogue of _composable_ packages for scientific computing. This enables `RelativisticDynamcis.jl` to easily leverage and interface with other scientific computing packages. In addition to providing a fast, modern package for strong field spin dynamics, `RelativisticDynamcis.jl` has two additional important features from the perspective of modern relativistic astrophysics. Firstly, it makes use of Julia's type-flexibility and is written in such a way so as to be able to support hardware accelerated low precision arithmetic and alternative rounding methods such as stochastic rounding. It is fully type flexible, to support arbitrary number formats for performance and analysis simultaneously. This means the model development is precision-agnostic, which allows us to address the common problems of dynamic range and critical precision loss often incurred from using low-precision number formats. This enables rapid prototyping and exploration towards the exploration of reduced precision numerical techniques in astrophysics, a technique common in other numerical fields such as climate modelling. Secondly, `RelativisticDynamcis.jl` is also written to be fully differentiable via automatic differentiation. This enables the package to be used for differentiable physics applications in astrophysics, for example gravitational waveform modelling and parameter estimation. The code is also written in a generally covariant way THIS IS ONLY MOSTLY TRUE (IF KERR - use analytical, otehrwise auto), with the translation to coordinate systems done implicitly. This makes further extensions to enable e.g. alternative spacetime metrics straightforward. Other useful features of `RelativisticDynamcis.jl` include the ability to specify initial conditions either in terms of conserved quantities (E,L,Q), which is useful when examining systems with very short orbital separations, typically used in the GW community or in terms of the familiar Keplerian elements which is more commonly used within the pulsar timing community. There is also a diagnostics tool for exxtracting gravitational waveforrs in the time domain via a numerical kuldge mehtod






particularly relevant for spinning systems like pulsars where the spin/rotation axis is very important 




The author and collaborators are have used the general methods and mathematics described in this package via legacy Fortran code, encoutnered multiple problems, and are now already using this package for multiple active
research projects



The author and collaborators have used the general methods and mathematics described in this package for multiple research projects, and this package represents an attempt to createa self documented well tested open source resource for 


Future extensios of this code would be to estend the dynamcis beyond second order in the multipole expansion, alternative spin conditions or alternative spacetime metrics. Since the code is written in a covariant way such extensions should be straightforward



higher order terms beyond the dipole 


progenitors of the extreme-mass-ratio-inspirals (EMRIs) and intermediate-mass-ratio-inspirals (IMRIs) gravitational wave events

via Fortran code, 


encoutnered multiple problems, and are now already using this package for multiple active
research projects



It solves the MPD equations numerically from a set of astrophysical initial conditions, allowing orbits to be specified in terms of the usual  In the cases where such codes exist, tehre is also a disconnect btween the orbital quantiteis such as XYZ and the conserved quanttie E,LQ







 primary focus is the dynamics - extensions include waveforms, unitary scaling for precision


is also written in a way 






whlst the research motivation in developing this package is in modeling relativistic pulsar systems, could also be extended generally to any systems with a sufficent mass ratio which need an accurate spin description. In the zero spo limit the solution reduces to the usual geodesic solution


unbound systems, parabolic etc.


Whilst the primary application is in modelling 




RelativisticDynamics.jl




 hC or Focombines 








# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References


# Scratch space




Its nature as a modern language also pr 


 . This enables  to easily interface with  and levarage 



modernness also enables easy extension, interfacing with other codes and import/levaaraging of other high performance libraries


moden also means dynamic, JIT types, package manager, simpel hardwar suport e.g. GPU

Enables interfacing with exixiting well tested packages much more readily than Fortran or C, enabling effective compatmentalisation of purpose and generality



Whilst the primary focus and research motivation is in application fast spinning MSP-BH systems, it reduces to the geodesic description in the zero spin case 

