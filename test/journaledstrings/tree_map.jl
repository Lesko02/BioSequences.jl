using DataStructures
include("journaledstrings.jl")


tree = AVLTree{Int64}()

for k in 1:2:20
    println(k)
    push!(tree, k)
end

println(tree)