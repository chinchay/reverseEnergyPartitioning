

def coef(N, V):
    """
    N: total number of particles
    V: volume of the system
    Perhaps you should change factorial(N-1) if not for a FCC. Read AppendixA
    """
    import math
    m = N - 1
    return math.factorial(m) / ( V ** m )
#

def belongs(x, xmin, xmax):
    return (xmin <= x) and (x <= xmax)
#

def productSqE(ε):
    import math
    from functools import reduce
    # http://book.pythontips.com/en/latest/map_filter.html
    product = reduce( (lambda x, y: x * y), ε )
    return 1 / math.sqrt( product )
#


def getForceMatrix():
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
    F = [ [f11, f12, f13],
          [f21, f22, f23],
          [f31, f32, f33]   ]
    return F
#

# function findEigen()
#     ε = eigvals( getForceMatrix() )
#     return ε
# end


def getRidZeros(ε):
    # harmonic oscillator has non-negative eigenvalues!
    # ε[k] != 0.0
    positives = []
    [positives.append(ε[k]) if ε[k] > 0 else None for k in range(len(ε)) ]
    return positives
#


def findEigenValues(M): # M is matrix array [ [], [], [], ...]
    import numpy as np
    from scipy import linalg as LA

    # convert array [1,2,3,] into numpy array: array([1,2,3,...])
    matrix = np.asarray(M)

    # For a complex Hermitian or real symmetric matrix: eigvalsh
    eigenValues = LA.eigvalsh(matrix)

    return eigenValues
#


# function to calculate the harmonic DOS
def harmonicDOS(V, dE):
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
    import math

    ε = findEigenValues( getForceMatrix()  )
    ε = getRidZeros(ε) # paper: "zeros do not contribute to DOS"
    D = len(ε) # D = 3N-3, but here it's not necessary to substract -3 since
                  # we already got rid of zeros.
    Dm = D / 2.0
    N  = (D + 3) / 3.0 # comes from solving D = 3N-3
    π = math.pi

    c  = coef(N, V)
    m1 = productSqE(ε) # not zeros!
    m2 = ( 2 * dE ) ** ( Dm  - 1 )
    m3 = 2 * ( π ** Dm ) / math.gamma(Dm)

    return c * m1 * m2 * m3
#
