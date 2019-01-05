include("functions.jl")

x = harmonicDOS(1, 1) #harmonicDOS(N, V, dE, D)


x = eigvals( getForceMatrix() )

x = getRidZeros(x)
println(x)


m2 = ( 2 * 1 ) ^ ( 1.5  - 1 )
m3 = 2 * ( π ^ 1.5 ) / gamma(1.5)

print(m2*m3)
