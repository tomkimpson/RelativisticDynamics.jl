# How to run SpeedyWeather.jl

The simplest way to run SpeedyWeather.jl with default parameters is

```julia
using SpeedyWeather
run_speedy()
```

Hooray, you have just simulated the Earth's atmosphere. Parameters, their meanings and
defaults can be found in `src/default_parameters.jl`, an incomplete list is provided below.
For example, if you want the simulation to run in double precision (`Float64`), at higher
resolution (`trunc`, the triangular spectral truncation), slow down the rotation of the Earth
(`rotation_earth` in ``s^{-1}``) and create some netCDF ouput, do

```julia
run_speedy(Float64,trunc=85,rotation_earth=1e-5,output=true)
```

If provided, the number format has to be the first argument, all other arguments are keyword arguments.
