var documenterSearchIndex = {"docs":
[{"location":"functions_index/#Index-of-functions","page":"Functions Index","title":"Index of functions","text":"","category":"section"},{"location":"functions_index/","page":"Functions Index","title":"Functions Index","text":"","category":"page"},{"location":"functions_index/","page":"Functions Index","title":"Functions Index","text":"Modules = [RelativisticDynamics]","category":"page"},{"location":"functions_index/#RelativisticDynamics.Constants","page":"Functions Index","title":"RelativisticDynamics.Constants","text":"C = Constants(P) A struct to hold all variables which are constant over the course of the integration. These are derived from the user-defined parameters \n\n\n\n\n\n","category":"type"},{"location":"functions_index/#RelativisticDynamics.Constants-Tuple{RelativisticDynamics.SystemParameters}","page":"Functions Index","title":"RelativisticDynamics.Constants","text":"Generator function for a Constants struct.\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.Model","page":"Functions Index","title":"RelativisticDynamics.Model","text":"\n\n\n\n","category":"type"},{"location":"functions_index/#RelativisticDynamics.PrognosticVariables","page":"Functions Index","title":"RelativisticDynamics.PrognosticVariables","text":"Struct holding the so-called 'prognostic' variables\n\n\n\n\n\n","category":"type"},{"location":"functions_index/#RelativisticDynamics.SystemParameters","page":"Functions Index","title":"RelativisticDynamics.SystemParameters","text":"P = Parameters(kwargs...) A struct to hold all model parameters that may be changed by the user. The struct uses keywords such that default values can be changed at creation. The default values of the keywords define the default model setup.\n\n\n\n\n\n","category":"type"},{"location":"functions_index/#RelativisticDynamics.ELQ-NTuple{5, Any}","page":"Functions Index","title":"RelativisticDynamics.ELQ","text":"Calculate the energy, angular momentum and Carter constant\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.Kretschmann_scalar-Tuple{Any, Any, Any}","page":"Functions Index","title":"RelativisticDynamics.Kretschmann_scalar","text":"Kretschman scalar for the Kerr metric\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.PlotTrajectory","page":"Functions Index","title":"RelativisticDynamics.PlotTrajectory","text":"Plot trajectory of a body. Assumes coordinates are Boyer Lindquist.  Plots in either 2D or 3D depending on specification of dimensions. Saves a low resolution PNG figure to disk in example_media/\n\n\n\n\n\n","category":"function"},{"location":"functions_index/#RelativisticDynamics.bounds_checks-Tuple{RelativisticDynamics.SystemParameters}","page":"Functions Index","title":"RelativisticDynamics.bounds_checks","text":"Look before you leap\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.christoffel-Tuple{Any, Any}","page":"Functions Index","title":"RelativisticDynamics.christoffel","text":"map = gridded(  alms::AbstractMatrix;\n                recompute_legendre::Bool=true,\n                grid::Type{<:AbstractGrid}=FullGaussianGrid)\n\nSpectral transform (spectral to grid space) from spherical coefficients alms to a newly allocated gridded field map. Based on the size of alms the grid type grid, the spatial resolution is retrieved based on the truncation defined for grid. SpectralTransform struct S is allocated to execute gridded(alms,S).\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.contravariant_metric-Tuple{Any, Any}","page":"Functions Index","title":"RelativisticDynamics.contravariant_metric","text":"Construct the NxN matrix of the contravariant metric. Metric components are defined via indvidual functions to allow for auto diff in unit tests\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.convert_to_covariant-Tuple{Any, Any}","page":"Functions Index","title":"RelativisticDynamics.convert_to_covariant","text":"Convert a vector from contravariant form to convariant form using the covariant metric \n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.covariant_metric-Tuple{Any, Any}","page":"Functions Index","title":"RelativisticDynamics.covariant_metric","text":"Construct the NxN matrix of the covariant metric. Metric components are defined via indvidual functions to allow for auto diff in unit tests\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.covariant_minkowski-Tuple{}","page":"Functions Index","title":"RelativisticDynamics.covariant_minkowski","text":"m = roundup_fft(n::Int;\n                small_primes::Vector{Int}=[2,3,5])\n\nReturns an integer m >= n with only small prime factors 2, 3, 5 (default, others can be specified with the keyword argument small_primes) to obtain an efficiently fourier-transformable number of longitudes, m = 2^i * 3^j * 5^k >= n, with i,j,k >=0.\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.delta-Tuple{Any, Any}","page":"Functions Index","title":"RelativisticDynamics.delta","text":"Δ = delta(r,a) The well-known delta function of the Kerr metric\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.initial_conditions-Tuple{RelativisticDynamics.Model}","page":"Functions Index","title":"RelativisticDynamics.initial_conditions","text":"Setup the initial conditions for the MPD orbital dynamics\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.mapping_d-Tuple{Any, Any, Any}","page":"Functions Index","title":"RelativisticDynamics.mapping_d","text":"d = mapping_d(r,a,zminus) Mapping function d used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.mapping_f-Tuple{Any, Any, Any}","page":"Functions Index","title":"RelativisticDynamics.mapping_f","text":"f = mapping_f(r,a,zminus) Mapping function f used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.mapping_g-Tuple{Any, Any}","page":"Functions Index","title":"RelativisticDynamics.mapping_g","text":"g = mapping_g(r,a) Mapping function g used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.mapping_h-Tuple{Any, Any, Any}","page":"Functions Index","title":"RelativisticDynamics.mapping_h","text":"h = mapping_h(r,a,zminus) Mapping function h used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.orbit-Union{Tuple{}, Tuple{Type{NF}}, Tuple{NF}} where NF<:AbstractFloat","page":"Functions Index","title":"RelativisticDynamics.orbit","text":"progn_vars = run_speedy(NF,kwargs...)\n\nRuns SpeedyWeather.jl with number format NF and any additional parameters in the keyword arguments kwargs.... Any unspecified parameters will use the default values as defined in src/parameters.jl.\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.riemann-Tuple{Any, Any}","page":"Functions Index","title":"RelativisticDynamics.riemann","text":"Riemann tensor components of the Kerr metric. First index is the contravariant, others are covariant   \n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.schwarzchild_covariant_riemann-Tuple{Any, Any}","page":"Functions Index","title":"RelativisticDynamics.schwarzchild_covariant_riemann","text":"Special case - the fully covariant components of the Riemann tensor for schwarzchild. From https://arxiv.org/pdf/0904.4184.pdf\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.sigma-Tuple{Any, Any, Any}","page":"Functions Index","title":"RelativisticDynamics.sigma","text":"Σ = sigma(r,θ,a) The well-known sigma function of the Kerr metric\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.spintensor-NTuple{4, Any}","page":"Functions Index","title":"RelativisticDynamics.spintensor","text":"Calculate the contravariant spin tensor. See e.g. https://mathworld.wolfram.com/PermutationTensor.html\n\n\n\n\n\n","category":"method"},{"location":"functions_index/#RelativisticDynamics.timestepping-Tuple{RelativisticDynamics.PrognosticVariables, RelativisticDynamics.Model}","page":"Functions Index","title":"RelativisticDynamics.timestepping","text":"The timesteppig Integration\n\n\n\n\n\n","category":"method"},{"location":"IC/#Initial-Conditions","page":"Initial Conditions","title":"Initial Conditions","text":"","category":"section"},{"location":"IC/","page":"Initial Conditions","title":"Initial Conditions","text":"The simplest way to run RelativisticDynamics.jl with default parameters is","category":"page"},{"location":"IC/","page":"Initial Conditions","title":"Initial Conditions","text":"using RelativisticDynamics\nsolution,model = orbit()","category":"page"},{"location":"IC/","page":"Initial Conditions","title":"Initial Conditions","text":"The orbit() funciton returns two objects. The first, solution holds the evolution of the position, momentum and spin vectors. The second, model, holds a copy of all the parameters and settings used to generate the solution (e.g. what was the BH spin?).","category":"page"},{"location":"IC/","page":"Initial Conditions","title":"Initial Conditions","text":"All default parameters can be found in src/default_parameters.jl. Passing a keyword argument to orbit() overrides the defaults e.g.","category":"page"},{"location":"IC/","page":"Initial Conditions","title":"Initial Conditions","text":"run_speedy(e=0.6,a=0.99)","category":"page"},{"location":"IC/","page":"Initial Conditions","title":"Initial Conditions","text":"would generate the solution for a system with an eccentricity = 0.6, around a BH with an extramal spin. ","category":"page"},{"location":"IC/","page":"Initial Conditions","title":"Initial Conditions","text":"If provided, the number format has to be the first argument, all other arguments are keyword arguments. e.g. ","category":"page"},{"location":"IC/","page":"Initial Conditions","title":"Initial Conditions","text":"run_speedy(Float32,e=0.6,a=0.99)","category":"page"},{"location":"IC/","page":"Initial Conditions","title":"Initial Conditions","text":"Please see notebooks/demo.ipynb for some worked examples using RelativisticDynamics.jl, including the application of autodiff methods.","category":"page"},{"location":"IC/","page":"Initial Conditions","title":"Initial Conditions","text":"<!– ","category":"page"},{"location":"IC/","page":"Initial Conditions","title":"Initial Conditions","text":"orbit","category":"page"},{"location":"IC/#RelativisticDynamics.orbit","page":"Initial Conditions","title":"RelativisticDynamics.orbit","text":"progn_vars = run_speedy(NF,kwargs...)\n\nRuns SpeedyWeather.jl with number format NF and any additional parameters in the keyword arguments kwargs.... Any unspecified parameters will use the default values as defined in src/parameters.jl.\n\n\n\n\n\n","category":"function"},{"location":"IC/","page":"Initial Conditions","title":"Initial Conditions","text":"–>","category":"page"},{"location":"how_to_run/#How-to-run-RelativisticDynamics.jl","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"","category":"section"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"The simplest way to run RelativisticDynamics.jl with default parameters is","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"using RelativisticDynamics\nsolution,model = orbit()","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"The orbit() funciton returns two objects. The first, solution holds the evolution of the position, momentum and spin vectors. The second, model, holds a copy of all the parameters and settings used to generate the solution (e.g. what was the BH spin?).","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"All default parameters can be found in src/default_parameters.jl. Passing a keyword argument to orbit() overrides the defaults e.g.","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"run_speedy(e=0.6,a=0.99)","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"would generate the solution for a system with an eccentricity = 0.6, around a BH with an extramal spin. ","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"If provided, the number format has to be the first argument, all other arguments are keyword arguments. e.g. ","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"run_speedy(Float32,e=0.6,a=0.99)","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"Please see notebooks/demo.ipynb for some worked examples using RelativisticDynamics.jl, including the application of autodiff methods.","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"<!– ","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"orbit","category":"page"},{"location":"how_to_run/#The-run_speedy-interface","page":"How to run RelativisticDynamics.jl","title":"The run_speedy interface","text":"","category":"section"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"orbit","category":"page"},{"location":"how_to_run/#The-run_speedy2-interface","page":"How to run RelativisticDynamics.jl","title":"The run_speedy2 interface","text":"","category":"section"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"bounds_checks","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"–>","category":"page"},{"location":"MPD/#Mathisson-Papetrou-Dixon-Equations","page":"MPD Equations","title":"Mathisson-Papetrou-Dixon Equations","text":"","category":"section"},{"location":"MPD/","page":"MPD Equations","title":"MPD Equations","text":"The simplest way to run RelativisticDynamics.jl with default parameters is","category":"page"},{"location":"MPD/","page":"MPD Equations","title":"MPD Equations","text":"using RelativisticDynamics\nsolution,model = orbit()","category":"page"},{"location":"MPD/","page":"MPD Equations","title":"MPD Equations","text":"The orbit() funciton returns two objects. The first, solution holds the evolution of the position, momentum and spin vectors. The second, model, holds a copy of all the parameters and settings used to generate the solution (e.g. what was the BH spin?).","category":"page"},{"location":"MPD/","page":"MPD Equations","title":"MPD Equations","text":"All default parameters can be found in src/default_parameters.jl. Passing a keyword argument to orbit() overrides the defaults e.g.","category":"page"},{"location":"MPD/","page":"MPD Equations","title":"MPD Equations","text":"run_speedy(e=0.6,a=0.99)","category":"page"},{"location":"MPD/","page":"MPD Equations","title":"MPD Equations","text":"would generate the solution for a system with an eccentricity = 0.6, around a BH with an extramal spin. ","category":"page"},{"location":"MPD/","page":"MPD Equations","title":"MPD Equations","text":"If provided, the number format has to be the first argument, all other arguments are keyword arguments. e.g. ","category":"page"},{"location":"MPD/","page":"MPD Equations","title":"MPD Equations","text":"run_speedy(Float32,e=0.6,a=0.99)","category":"page"},{"location":"MPD/","page":"MPD Equations","title":"MPD Equations","text":"Please see notebooks/demo.ipynb for some worked examples using RelativisticDynamics.jl, including the application of autodiff methods.","category":"page"},{"location":"MPD/","page":"MPD Equations","title":"MPD Equations","text":"orbit","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = RelativisticDynamics","category":"page"},{"location":"#RelativisticDynamics.jl-documentation","page":"Home","title":"RelativisticDynamics.jl documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation for RelativisticDynamics.jl, a relativistic orbital dynamics model for simulating spinning objects in curved spacetime.","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"RelativisticDynamics.jl is a numerical model for determining the spin-orbital dynamics of a spinning body, such as a pulsar, on a background Kerr spacetime. The code solves a set of ODEs numerically. These equations are based on the original works of Mathisson 1937, Papapetrou 1951 and Dixon 1964. Consequently these equations are known as the MPD equations. More recent works can be found in Mashoon & Singh 2006, Singh 2005, Singh, Wu & Sarty 2014 and Li,Wu & Singh 2019. Additional interesting discussion on the motion of extended bodies in GR can be found in Costa & Natário, 2015","category":"page"},{"location":"#Manual-Outline","page":"Home","title":"Manual Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Please see the following pages of the documentation for more details    ","category":"page"},{"location":"","page":"Home","title":"Home","text":"How to run RelativisticDynamics.jl\nInitial conditions\nMPD Equations\nIndex of functions","category":"page"},{"location":"#Scope","page":"Home","title":"Scope","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The focus of SpeedyWeather.jl is to develop a global atmospheric model of intermediate complexity, that can run at various levels of precision (16, 32 and 64-bit) on different architectures (x86 and ARM, currently planned, GPUs probably in the future). Additionally, the model is written in an entirely number format-flexible way, such that any custom number format can be used and Julia will compile to the format automatically.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SpeedyWeather.jl is registered in the Julia Registry. Open Julia's package manager from the REPL with ] and add the github repository to install SpeedyWeather.jl and all dependencies","category":"page"},{"location":"","page":"Home","title":"Home","text":"(@v1.7) pkg> add SpeedyWeather","category":"page"},{"location":"","page":"Home","title":"Home","text":"which will automatically install the latest release. However, you may want to install directly from the main branch with","category":"page"},{"location":"","page":"Home","title":"Home","text":"(@v1.7) pkg> add https://github.com/milankl/SpeedyWeather.jl#main","category":"page"},{"location":"","page":"Home","title":"Home","text":"other branches than #main can be installed by adding #branch_name instead.","category":"page"},{"location":"#Acknowledgements","page":"Home","title":"Acknowledgements","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This work solving the MPD equations was originally motivated through interesting discussions with Kinwah Wu. The port to a modern, precision-flexible model in Julia was heavily inspired by Milan Klöwer. Huge thanks to both.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Contributions are always welcome - just open a pull request","category":"page"}]
}
