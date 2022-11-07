var documenterSearchIndex = {"docs":
[{"location":"how_to_run/#How-to-run-SpeedyWeather.jl","page":"How to run RelativisticDynamics.jl","title":"How to run SpeedyWeather.jl","text":"","category":"section"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"The simplest way to run SpeedyWeather.jl with default parameters is","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"using SpeedyWeather\nrun_speedy()","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"Hooray, you have just simulated the Earth's atmosphere. Parameters, their meanings and defaults can be found in src/default_parameters.jl, an incomplete list is provided below. For example, if you want the simulation to run in double precision (Float64), at higher resolution (trunc, the triangular spectral truncation), slow down the rotation of the Earth (rotation_earth in s^-1) and create some netCDF ouput, do","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"run_speedy(Float64,trunc=85,rotation_earth=1e-5,output=true)","category":"page"},{"location":"how_to_run/","page":"How to run RelativisticDynamics.jl","title":"How to run RelativisticDynamics.jl","text":"If provided, the number format has to be the first argument, all other arguments are keyword arguments.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = RelativisticDynamics","category":"page"},{"location":"#RelativisticDynamics.jl-documentation","page":"Home","title":"RelativisticDynamics.jl documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation for RelativisticDynamics.jl a relativistic orbital dynamics model for simulating spinning objects in curved spacetime.","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"RelativisticDynamics.jl is a numerical model for determining the spin-orbital dynamics of a spinning body, such as a pulsar, on a background Kerr spacetime. The code solves a set of ODEs numerically. These equations are based on the original works of Mathisson 1937, Papapetrou 1951 and Dixon 1964. Consequently these equations are known as the MPD equations. More recent works can be found in Mashoon & Singh 2006, Singh 2005, Singh, Wu & Sarty 2014 and Li,Wu & Singh 2019. Additional interesting discussion on the motion of extended bodies in GR can be found in Costa & Natário, 2015","category":"page"},{"location":"#Manual-Outline","page":"Home","title":"Manual Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Please see the following pages of the documentation for more details    ","category":"page"},{"location":"","page":"Home","title":"Home","text":"How to run RelativisticDynamics.jl\nMPD Equations\nIndex of functions","category":"page"},{"location":"#Scope","page":"Home","title":"Scope","text":"","category":"section"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#Acknowledgements","page":"Home","title":"Acknowledgements","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This work solving the MPD equations was originally motivated through interesting discussions with Kinwah Wu. The port to a modern, precision-flexible model in Julia was heavily inspired by Milan Klöwer. Huge thanks to both.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Contributions are always welcome - just open a pull request","category":"page"},{"location":"#THIS-HEADING","page":"Home","title":"THIS HEADING","text":"","category":"section"},{"location":"#THAT-HEADING","page":"Home","title":"THAT HEADING","text":"","category":"section"},{"location":"#SOMETHING-ELSE","page":"Home","title":"SOMETHING ELSE","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [RelativisticDynamics]","category":"page"},{"location":"#RelativisticDynamics.Constants","page":"Home","title":"RelativisticDynamics.Constants","text":"C = Constants(P) A struct to hold all variables which are constant over the course of the integration. These are derived from the user-defined parameters \n\n\n\n\n\n","category":"type"},{"location":"#RelativisticDynamics.Constants-Tuple{RelativisticDynamics.SystemParameters}","page":"Home","title":"RelativisticDynamics.Constants","text":"Generator function for a Constants struct.\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.Model","page":"Home","title":"RelativisticDynamics.Model","text":"\n\n\n\n","category":"type"},{"location":"#RelativisticDynamics.PrognosticVariables","page":"Home","title":"RelativisticDynamics.PrognosticVariables","text":"Struct holding the so-called 'prognostic' variables\n\n\n\n\n\n","category":"type"},{"location":"#RelativisticDynamics.SystemParameters","page":"Home","title":"RelativisticDynamics.SystemParameters","text":"P = Parameters(kwargs...) A struct to hold all model parameters that may be changed by the user. The struct uses keywords such that default values can be changed at creation. The default values of the keywords define the default model setup.\n\n\n\n\n\n","category":"type"},{"location":"#RelativisticDynamics.ELQ-NTuple{6, Any}","page":"Home","title":"RelativisticDynamics.ELQ","text":"Calculate the energy, angular momentum and Carter constant for the different types of system\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.Kretschmann_scalar-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.Kretschmann_scalar","text":"Kretschman scalar for the Kerr metric\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.MPD_initial_conditions-Tuple{RelativisticDynamics.Model}","page":"Home","title":"RelativisticDynamics.MPD_initial_conditions","text":"Setup the initial conditions for the MPD orbital dynamics\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.PlotTrajectory","page":"Home","title":"RelativisticDynamics.PlotTrajectory","text":"Plot trajectory of a body. Assumes coordinates are Boyer Lindquist.  Plots in either 2D or 3D depending on specification of dimensions. Saves a low resolution PNG figure to disk in example_media/\n\n\n\n\n\n","category":"function"},{"location":"#RelativisticDynamics.bounds_checks-Tuple{RelativisticDynamics.SystemParameters}","page":"Home","title":"RelativisticDynamics.bounds_checks","text":"Look before you leap\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.christoffel-Tuple{Any, Any}","page":"Home","title":"RelativisticDynamics.christoffel","text":"map = gridded(  alms::AbstractMatrix;\n                recompute_legendre::Bool=true,\n                grid::Type{<:AbstractGrid}=FullGaussianGrid)\n\nSpectral transform (spectral to grid space) from spherical coefficients alms to a newly allocated gridded field map. Based on the size of alms the grid type grid, the spatial resolution is retrieved based on the truncation defined for grid. SpectralTransform struct S is allocated to execute gridded(alms,S).\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.contravariant_metric-Tuple{Any, Any}","page":"Home","title":"RelativisticDynamics.contravariant_metric","text":"Construct the NxN matrix of the contravariant metric. Metric components are defined via indvidual functions to allow for auto diff in unit tests\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.convert_to_covariant-Tuple{Any, Any}","page":"Home","title":"RelativisticDynamics.convert_to_covariant","text":"Convert a vector from contravariant form to convariant form using the covariant metric \n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.covariant_minkowski-Tuple{}","page":"Home","title":"RelativisticDynamics.covariant_minkowski","text":"Construct the NxN matrix of the covariant Minkowski metric.\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.delta-Tuple{Any, Any}","page":"Home","title":"RelativisticDynamics.delta","text":"Δ = delta(r,a) The well-known delta function of the Kerr metric\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.mapping_d-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.mapping_d","text":"d = mapping_d(r,a,zminus) Mapping function d used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.mapping_f-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.mapping_f","text":"f = mapping_f(r,a,zminus) Mapping function f used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.mapping_g-Tuple{Any, Any}","page":"Home","title":"RelativisticDynamics.mapping_g","text":"g = mapping_g(r,a) Mapping function g used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.mapping_h-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.mapping_h","text":"h = mapping_h(r,a,zminus) Mapping function h used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.orbit-Union{Tuple{}, Tuple{Type{NF}}, Tuple{NF}} where NF<:AbstractFloat","page":"Home","title":"RelativisticDynamics.orbit","text":"The main output function from this package.\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.riemann-Tuple{Any, Any}","page":"Home","title":"RelativisticDynamics.riemann","text":"Riemann tensor components of the Kerr metric. First index is the contravariant, others are covariant   \n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.schwarzchild_covariant_riemann-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.schwarzchild_covariant_riemann","text":"Special case - the fully covariant components of the Riemann tensor for schwarzchild. From https://arxiv.org/pdf/0904.4184.pdf\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.sigma-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.sigma","text":"Σ = sigma(r,θ,a) The well-known sigma function of the Kerr metric\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.spintensor-NTuple{4, Any}","page":"Home","title":"RelativisticDynamics.spintensor","text":"Calculate the contravariant spin tensor. See e.g. https://mathworld.wolfram.com/PermutationTensor.html\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.timestepping-Tuple{RelativisticDynamics.PrognosticVariables, RelativisticDynamics.Model}","page":"Home","title":"RelativisticDynamics.timestepping","text":"The timesteppig Integration\n\n\n\n\n\n","category":"method"}]
}
