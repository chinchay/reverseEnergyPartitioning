import Pkg
Pkg.add("SpecialFunctions")
using SpecialFunctions # for gamma() function
using LinearAlgebra # for eigvals()


function coef(N, V)
    """
    N: total number of particles
    V: volume of the system
    Perhaps you should change factorial(N-1) if not for a FCC. Read AppendixA
    """
    m = N - 1
    return factorial(m) / ( V ^ m )
end

function belongs(x, xmin, xmax)
    return (xmin <= x) & (x <= xmax)
end

function productSqE(ε)
    δ = 1.0E-5
    p = 1.0
    for k in 1:1:length(ε)
        @assert ε[k] >= 0.0 # harmonic oscillator has non-negative eigenvalues!
        if !belongs(ε[k], 0.0, δ)  # ε[k] != 0.0
            p = p / sqrt( ε[k] )  # ( ε[k] ^ ( -0.5 ) )
        end
    end
    return p
end

function getForceMatrix()
    """
    This function get the Hessian
    """
    f11 = 1.0
    f12 = 0.0
    f13 = 0.0

    f21 = 0.0
    f22 = 1.0
    f23 = 0.0

    f31 = 0.0
    f32 = 0.0
    f33 = 1.0

    # F = [ 0.0 1.0;
    #       1.0 0.0 ]

    #
    F = [ f11 f12 f13;
          f21 f22 f23;
          f31 f32 f33   ]
    return F

end

# function findEigen()
#     ε = eigvals( getForceMatrix() )
#     return ε
# end

function getRidZeros(ε)
    δ = 1.0E-5
    v = zeros(0) # This helps to initialize with type float, but length ZERO!
    for k in 1:1:length(ε)
        @assert ε[k] >= 0.0 # harmonic oscillator has non-negative eigenvalues!
        if !belongs(ε[k], 0.0, δ)  # ε[k] != 0.0
            push!(v, ε[k])
        end
    end
    return v
end

# function to calculate the harmonic DOS
function harmonicDOS(V, dE)
    """
    Calculate the harmonic DOS
    Here one particle is taken as center of reference, and the other particles
    are explained with respect to it. For two particles joined by a string in a
    1-dimensional case:
    V(x1,x2) = k12 * (1/2) * (x1 - x2)^2
    V(x1,x2) = k12 * (1/2) * (x1^2 - 2*x1*x2 + x2^2)
    then the Hessian (force matrix) is:
    F11 = k12 * (1/2) * (2)
    F12 = k12 * (1/2) * (-2)
    F22 = k12 * (1/2) * (2)
    F = k12 [  1  -1;
              -1   1  ]
    with eigenvalues {0, 2}. We would obtain eigenvalue=2 getting rid of the
    "translational normal modes" by describing x1 with respect to x2:
    V(r21) = k12 * r12^2
    The expression V(x1,x2) is a particular case of the rotated parabole
    equation:
    Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
    rotated by an angle θ that tan(2θ) = B / (A-C)
    A particular case: V(x1,x2) = A*x1^2 - 2*sqrt(A*C) + C*x2^2
                                = (A*x1 - B*x2)^2
    which can represent the interaction of two particles joined by a string.

    """
    ε = eigvals( getForceMatrix() )
    ε = getRidZeros(ε) # paper: "zeros do not contribute to DOS"
    D = length(ε) # D = 3N-3, but here it's not necessary to substract -3 since
                  # we already got rid of zeros.
    Dm = D / 2.0
    N  = (D + 3) / 3.0 # comes from solving D = 3N-3

    c  = coef(N, V)
    m1 = productSqE(ε)
    m2 = ( 2 * dE ) ^ ( Dm  - 1 )
    m3 = 2 * ( π ^ Dm ) / gamma(Dm)

    return c * m1 * m2 * m3

end
