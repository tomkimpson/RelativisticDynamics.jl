var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = RelativisticDynamics","category":"page"},{"location":"#RelativisticDynamics","page":"Home","title":"RelativisticDynamics","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for RelativisticDynamics.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation for RelativisticDynamics.jl a relativistic orbital dynamics model for simulating spinning objects in curved spacetime.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [RelativisticDynamics]","category":"page"},{"location":"#RelativisticDynamics.Constants","page":"Home","title":"RelativisticDynamics.Constants","text":"C = Constants(P)\n\nA struct to hold all variables which are constant over the course of the integration.\n\n\n\n\n\n","category":"type"},{"location":"#RelativisticDynamics.Constants-Tuple{SystemParameters}","page":"Home","title":"RelativisticDynamics.Constants","text":"Generator function for a Constants struct.\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.Model","page":"Home","title":"RelativisticDynamics.Model","text":"M = BarotropicModel(::Parameters,\n                    ::Constants,\n                    ::GeoSpectral,\n                    ::HorizontalDiffusion)\n\nThe BarotropicModel struct holds all other structs that contain precalculated constants, whether scalars or arrays that do not change throughout model integration. In contrast to ShallowWaterModel or PrimitiveEquationModel it does not contain a Boundaries struct as not needed.\n\n\n\n\n\n","category":"type"},{"location":"#RelativisticDynamics.PrognosticVariables","page":"Home","title":"RelativisticDynamics.PrognosticVariables","text":"Struct holding the so-called 'prognostic' variables\n\n\n\n\n\n","category":"type"},{"location":"#RelativisticDynamics.SystemParameters","page":"Home","title":"RelativisticDynamics.SystemParameters","text":"P = Parameters(kwargs...)\n\nA struct to hold all model parameters that may be changed by the user. The struct uses keywords such that default values can be changed at creation. The default values of the keywords define the default model setup.\n\n\n\n\n\n","category":"type"},{"location":"#RelativisticDynamics.delta-Tuple{Any, Any}","page":"Home","title":"RelativisticDynamics.delta","text":"Δ = delta(r,a) The well-known function of the Kerr metric\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.initial_conditions-Tuple{RelativisticDynamics.Model}","page":"Home","title":"RelativisticDynamics.initial_conditions","text":"Initialize a PrognosticVariables struct for an atmosphere at rest. No winds, hence zero vorticity and divergence, but temperature, pressure and humidity are initialised \n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.mapping_d-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.mapping_d","text":"d = mapping_d(r,a,zminus) Mapping function d used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.mapping_f-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.mapping_f","text":"f = mapping_f(r,a,zminus) Mapping function f used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.mapping_g-Tuple{Any, Any}","page":"Home","title":"RelativisticDynamics.mapping_g","text":"g = mapping_g(r,a) Mapping function g used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.mapping_h-Tuple{Any, Any, Any}","page":"Home","title":"RelativisticDynamics.mapping_h","text":"h = mapping_h(r,a,zminus) Mapping function h used when converting from Keplerian orbital parameters to constants of motion\n\n\n\n\n\n","category":"method"},{"location":"#RelativisticDynamics.orbit-Union{Tuple{}, Tuple{Type{NF}}, Tuple{NF}} where NF<:AbstractFloat","page":"Home","title":"RelativisticDynamics.orbit","text":"Some docstring for the run_program function\n\n\n\n\n\n","category":"method"}]
}
