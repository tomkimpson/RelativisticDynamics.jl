using LaTeXStrings
using RecipesBase


"""
Plotting recipe for use with Plots.jl
"""
@recipe function f(sol::ODESolution,a::Float64; upsample=10,vars_to_plot=[:x,:y])

    #Default plotting settings
    size      -->  (600, 600)
    linecolor --> 1
    legend    --> :none

    #Up-sample the solution according to the upsample factor 
    T = range(first(sol.t),last(sol.t),length=length(sol.t)*upsample)
    p = sol(T)

    # Extract relevant data from the interpolated solution 
    r = p[2,:]
    θ = p[3,:] 
    ϕ = p[4,:]

    # Boyer lindquist to Cartesian 
    w = sqrt.(r.^2 .+ a^2) 
    x = w .* sin.(θ) .* cos.(ϕ)
    y = w .* sin.(θ) .* sin.(ϕ)
    z = r .* cos.(θ)


    #Define a basic dict that accepts a :symbol key and returns an upsampled solution for that variable + the axis label
    #This is a bit verbose - there must be a more concise way, but for our smll number of vars this is OK for now.
    vars_dict = Dict(
                    :x => [x,L"x (r_h)"],
                    :y => [y,L"y (r_h)"],
                    :z => [z,L"z (r_h)"],

                    :t => [T,L"t"],
                    :r => [r,L"r"],
                    :θ => [θ,L"θ"],
                    :ϕ => [ϕ,L"ϕ"],

                    :pt => [p[5,:],L"p^t"],
                    :pr => [p[6,:],L"p^r"],
                    :pθ => [p[7,:],L"p^{θ}"],
                    :pϕ => [p[8,:],L"p^{ϕ}"],
                    
                    :st => [p[9,:],L"s^t"],
                    :sr => [p[10,:],L"s^r"],
                    :sθ => [p[11,:],L"s^{θ}"],
                    :sϕ => [p[12,:],L"s^{ϕ}"],
    
    )

    #Axis labels
    xlabel    --> vars_dict[vars_to_plot[1]][2]
    ylabel    --> vars_dict[vars_to_plot[2]][2]

    #data to plot
    vars_dict[vars_to_plot[1]][1],vars_dict[vars_to_plot[2]][1]
end 








