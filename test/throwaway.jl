using TensorOperations



pvector = zeros(Float64,4)
Riemann = zeros(Float64,4,4,4,4)
Stensor = zeros(Float64,4,4)
Riemann
m0 = 1.0



@tensor begin
    dx[α] := pvector[α]+0.50*Stensor[α,β]*Riemann[β,γ,μ,ν]*pvector[γ]*Stensor[μ,ν]/(m0^2)
end