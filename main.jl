include("functions.jl")
println(harmonicDOS(1, 1))


include("potential.jl")

atoms()

atom1 = Atom("Co", 0.0, 0.0)

sentence = "This is a string"
sArray = split(sentence)
longest = reduce( (x,y) -> length(x) > length(y) ? x : y, sArray )
print(longest)

v = [1000,1,2,3,-1,4,100,99,-99,101]
m = reduce( (x,y) -> x > y ? x : y, v )


reduce(push, 1, [2; 3; 4])

Ω
α


struct Individual{T}
  chromosome::Vector{T}

  #Inner constructor
  Individual{T}(ngenes::UInt) where T<:Number =  new{T}(Vector{T}(ngenes))
end


struct BinaryIndividual <: Individual

  #Inner constructor
  BinaryIndividual(ngenes::UInt) = new{Bool}(ngenes)

end
