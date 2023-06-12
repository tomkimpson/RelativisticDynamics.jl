

using Zygote


function some_array(coords,a) # `coords` is a vector length 4, `a` is a scalar 

    xs = zeros(typeof(a),4,4)
    g = Zygote.bufferfrom(xs) #a 4x4 buffered Zygote object

    t,r,θ,ϕ =  coords[1],coords[2],coords[3],coords[4] #unpack the coordinates

    g[1,1] = 2*r #assign array element

    return copy(g)


end 

#Coordinates [t,r,θ,ϕ] and parameter `a` 
t = 0.0
r = 10.0
θ = π/2.0 
ϕ = 0.0 
a = 0.1
coords = [t,r,θ,ϕ]

#Create an array 
some_array(coords,a)

#Get jacobian of array. Works
array_jacobian = jacobian(x -> some_array(x,a), coords)[1]
array_jacobian = reshape(array_jacobian,(4,4,4)) #array_jacobian[:, :, i] is the derivative of some_array w.r.t. coords[i]

#Get hessian of array. Does not work?
array_hessian = hessian(x -> some_array(x,a), coords)[1]