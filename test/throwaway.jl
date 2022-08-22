using TensorOperations

p = [-1.4855661806711106e9, -6.422023154178244e8, -4.292067380342278e6, -1.150055988422677e6]
#p = [-1, -6, 6, 5]

g = zeros(Float64,4,4)
# g[1,1] = -0.9933333333333333
# g[2,2] = 1.0067114093959733
# g[3,3] = 90000.0
# g[4,4] = 90000.0

g[1,1] = -0.9933333333333
g[2,2] = 1.0067114093959
g[3,3] = 90000.0
g[4,4] = 90000.0




display(g)

@tensor begin 
    Vsq = g[μ,ν]*p[μ]*p[ν]
end 

println("TensorOperations result:")
println(Vsq)


println("Manual result")
man_res = g[1,1]*p[1]*p[1] + g[2,2]*p[2]*p[2] + g[3,3]*p[3]*p[3] + g[4,4]*p[4]*p[4]
println(man_res)