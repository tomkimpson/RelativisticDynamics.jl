---
title: 'RelativisticDynamics.jl: Relativisitc Spin-Orbital Dynamics in Julia'
tags:
  - Julia
  - astronomy
  - dynamics
  - gravity
  - general relativity
authors:
  - name: Tom Kimpson 1
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: 1, 2 # (Multiple affiliations must be quoted)
  - name: Tom Kimpson 2
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
affiliations:
 - name: School of Physics, University of Melbourne, Parkville, VIC 3010, Australia
   index: 1
 - name: OzGrav, University of Melbourne, Parkville, VIC 3010, Australia
   index: 2
date: 14 October 2022
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary



Relativistic binaries composed of a millisecond pulsar (MSP) orbiting a much more massive BH ($\gtrsim 10^3 M_{\odot}$) are exceptional probes for investigating key questions of fundamental physics and astrophysics. 


Such systems are natural sources of gravitational waves in the mHz regime, expected to be detectable by LISA. The associated radio emission from the companion pulsar raises the possibility of an EM counterpart enabling high precision multimessenger measurements to be made.  



The description of the orbital dynamics of these systems, and the influence on the resultant EM and GW observed signal, is non-trivial. Whilst models commonly treat the motion via a PN or geodesic description, such an approach neglects the non-linear influence of the pulsar spin on the underlying spacetime metric.

A proper treatment of the spin-orbital dynamics can be derived from the conservation of the energy-momentum tensor
\begin{equation}\label{eq:conservation}
\hat T^{\mu \nu}_{;\nu} = 0
\end{equation}
which when expanded into a set of infinite multipole moments leads to a description of the momentum vector $p^{\mu}$ and the spin tensor $s^{\mu \nu}$ 



\begin{equation}\label{eq:mpd1}
\hat \frac{Dp^{\mu}}{d \lambda} = -\frac{1}{2}R^{\mu}_{\nu \alpha \beta} u^{\nu} s^{\alpha \beta}
\end{equation}

\begin{equation}\label{eq:mpd2}
\hat \frac{Ds^{\mu \nu}}{d \lambda} =p^{\mu}u^{\nu} - p^{\nu}u^{\mu}
\end{equation}


Numerical



These equations togetehr describe the spin-orbital evolution in a fullly conssiten way that is applicable to trong field regimes, and are knownas the MPD equations



whlst the research motivation in developing this package is in modeling relativistic pulsar systems, could also be extended generally to any systems with a sufficent mass ratio which need an accurate spin description. In the zero spo limit the solution reduces to the usual geodesic solution





The author and collaborators are have used the general methods and mathematics described in this package via legacy Fortran code, encoutnered multiple problems, and are now already using this package for multiple active
research projects

Whilst the primary application is in modelling 

RelativisticDynamics.jl is an open source code for numerical evolution of spinning relativistic systems 

# Statement of need

`RelativisticDynamics.jl` is Julia package for relativistic spin orbital dynamics in the gravitational strong field. Julia is a modern language that solves the "two language problem" CITE, enabling fast typing with petaflop performance, comparable with better known numerical lanuages such as C or Fortran.




Existing codes for modelling spinning bodies in the gravitational strong field are generally lacking; at the low precisin end newtonian dynakics with PN corrections, o geodesic orbital dynamics whilst a full numerical relativity evolution is unnecesaryly computaitonlly and untractable for these systems which are observed over a large numberof orbital cycles intensive for these class of extreme mass ratio problems where the smaller body can be considered as a test mass


In the cases where such codes exist, tehre is also a disconnect btween the orbital quantiteis such as XYZ and the conserved quanttie E,LQ


modernness also enables easy extension, interfacing with other codes and import/levaaraging of other high performance libraries


moden also means dynamic, JIT types, package manager, simpel hardwar suport e.g. GPU

Enables interfacing with exixiting well tested packages much more readily than Fortran or C, enabling effective compatmentalisation of purpose and generality



In addition to providing a fast, modern package for strong field spin dynamics, Relativistic dynamics has two additional importnat features from the perspective of modern relativistic astrophysics. Firstly it makes use of ulias type-flexibility and is written in such a way so as to be able ti support hardware accenerlated low preision artihmetic and alternaive rounding methids such as stchastic rounding (cite MILAN). It is fully type flexible, to support arbitrary number formats for performance and analysis simultaneously. This means the model development is precision-agnostic, which allows us to address the common problems of dynamic range and critical precision loss often incurred from using low-precision number formats. This enables rapid prototyping and exploration towards the exploration of reduced precision numerical techniques in astrophysics, where they are already commonly used in e.g. numerical climate simulations.



Relativistic dynamcis is also written to be fully differentialble via autoamtic differentiation. This raises the possobility of using such a package for differentiable physics applicatins in astrophysics e.g. gravitational waveform modelling and parameter estimation. 


is also written in a way 





RelativisticDynamics.jl




 hC or Focombines 



Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations

and refer to \autoref{eq:fourier} from text.

# Key features?




# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

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