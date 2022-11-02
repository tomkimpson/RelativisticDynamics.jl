var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = RelativisticDynamics","category":"page"},{"location":"#RelativisticDynamics","page":"Home","title":"RelativisticDynamics","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for RelativisticDynamics.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation for RelativisticDynamics.jl a relativistic orbital dynamics model for simulating spinning objects in curved spacetime.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [RelativisticDynamics]","category":"page"},{"location":"#RelativisticDynamics.Constants","page":"Home","title":"RelativisticDynamics.Constants","text":"C = Constants(P) A struct to hold all variables which are constant over the course of the integration. These are derived from the user-defined parameters \n\n\n\n\n\n","category":"type"},{"location":"#RelativisticDynamics.Constants-Tuple{SystemParameters}","page":"Home","title":"RelativisticDynamics.Constants","text":"Generator function for a Constants struct.\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.Model","page":"Home","title":"RelativisticDynamics.Model","text":"\n\n\n\n","category":"type"},{"location":"#RelativisticDynamics.PrognosticVariables","page":"Home","title":"RelativisticDynamics.PrognosticVariables","text":"Struct holding the so-called 'prognostic' variables\n\n\n\n\n\n","category":"type"},{"location":"#RelativisticDynamics.SystemParameters","page":"Home","title":"RelativisticDynamics.SystemParameters","text":"P = Parameters(kwargs...) A struct to hold all model parameters that may be changed by the user. The struct uses keywords such that default values can be changed at creation. The default values of the keywords define the default model setup.\n\n\n\n\n\n","category":"type"},{"location":"#RelativisticDynamics.ELQ-NTuple{6, Any}","page":"Home","title":"RelativisticDynamics.ELQ","text":"Calculate the energy, angular momentum and Carter constant for the different types of system\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.Kretschmann_scalar-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.Kretschmann_scalar","text":"Kretschman scalar for the Kerr metric\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.MPD_initial_conditions-Tuple{RelativisticDynamics.Model}","page":"Home","title":"RelativisticDynamics.MPD_initial_conditions","text":"Setup the initial conditions for the MPD orbital dynamics\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.PlotTrajectory","page":"Home","title":"RelativisticDynamics.PlotTrajectory","text":"Plot trajectory of a body. Assumes coordinates are Boyer Lindquist.  Plots in either 2D or 3D depending on specification of dimensions. Saves a low resolution PNG figure to disk in example_media/\n\n\n\n\n\n","category":"function"},{"location":"#RelativisticDynamics.bounds_checks-Tuple{SystemParameters}","page":"Home","title":"RelativisticDynamics.bounds_checks","text":"Look before you leap\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.christoffel-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.christoffel","text":"Christoffel symbols of the second kind \n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.contravariant_metric-Tuple{Any, Any}","page":"Home","title":"RelativisticDynamics.contravariant_metric","text":"Construct the NxN matrix of the contravariant metric. Metric components are defined via indvidual functions to allow for auto diff in unit tests\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.convert_to_covariant-Tuple{Any, Any}","page":"Home","title":"RelativisticDynamics.convert_to_covariant","text":"Convert a vector from contravariant form to convariant form using the covariant metric \n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.covariant_metric-Tuple{Any, Any}","page":"Home","title":"RelativisticDynamics.covariant_metric","text":"Construct the NxN matrix of the covariant metric. Metric components are defined via indvidual functions to allow for auto diff in unit tests\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.covariant_minkowski-Tuple{}","page":"Home","title":"RelativisticDynamics.covariant_minkowski","text":"Construct the NxN matrix of the covariant Minkowski metric.\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.delta-Tuple{Any, Any}","page":"Home","title":"RelativisticDynamics.delta","text":"Δ = delta(r,a) The well-known delta function of the Kerr metric\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.mapping_d-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.mapping_d","text":"d = mapping_d(r,a,zminus) Mapping function d used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.mapping_f-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.mapping_f","text":"f = mapping_f(r,a,zminus) Mapping function f used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.mapping_g-Tuple{Any, Any}","page":"Home","title":"RelativisticDynamics.mapping_g","text":"g = mapping_g(r,a) Mapping function g used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.mapping_h-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.mapping_h","text":"h = mapping_h(r,a,zminus) Mapping function h used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.orbit-Union{Tuple{}, Tuple{Type{NF}}, Tuple{NF}} where NF<:AbstractFloat","page":"Home","title":"RelativisticDynamics.orbit","text":"The main output function from this package.\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.riemann-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.riemann","text":"Riemann tensor components of the Kerr metric. First index is the contravariant, others are covariant   \n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.schwarzchild_covariant_riemann-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.schwarzchild_covariant_riemann","text":"Special case - the fully covariant components of the Riemann tensor for schwarzchild. From https://arxiv.org/pdf/0904.4184.pdf\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.sigma-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.sigma","text":"Σ = sigma(r,θ,a) The well-known sigma function of the Kerr metric\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.spintensor-NTuple{4, Any}","page":"Home","title":"RelativisticDynamics.spintensor","text":"Calculate the contravariant spin tensor. See e.g. https://mathworld.wolfram.com/PermutationTensor.html\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.timestepping-Tuple{PrognosticVariables, RelativisticDynamics.Model}","page":"Home","title":"RelativisticDynamics.timestepping","text":"The timesteppig Integration\n\n\n\n\n\n","category":"method"}]
}
