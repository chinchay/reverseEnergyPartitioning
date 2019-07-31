

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

def productSqE(eig):
    import math
    from functools import reduce
    # http://book.pythontips.com/en/latest/map_filter.html
    product = reduce( (lambda x, y: x * y), eig )
    return 1 / math.sqrt( product )
#


def getForceMatrix():
    """
    This function get the Hessian
    """
    import numpy as np
    nAtoms = 2
    F = np.diag(np.ones(3 * nAtoms), 0) # each atom as x,y,z degree of freedom

    return F #returns a numpy array
#




def getRidZeros(eig):
    # harmonic oscillator has non-negative eigenvalues!
    # eig[k] != 0.0
    positives = []
    [positives.append(eig[k]) if eig[k] > 0 else None for k in range(len(eig)) ]
    return positives
#


# def findEigenValues(M): # M is numpy matrix array [ [], [], [], ...]
#     # to convert array to numpy array, use np.asarray( [array] )
#     # M is already a numpy array.
#     import numpy as np
#     from scipy import linalg as LA
#
#     # For a complex Hermitian or real symmetric matrix: eigvalsh
#     eigenValues = LA.eigvalsh(M)
#
#     return eigenValues
# #


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
    rotated by an angle theta that tan(2theta) = B / (A-C)
    A particular case: V(x1,x2) = A*x1^2 - 2*sqrt(A*C) + C*x2^2
                                = (A*x1 - B*x2)^2
    which can represent the interaction of two particles joined by a string.

    """
    import math
    from scipy import linalg as LA

    # For a complex Hermitian or real symmetric matrix: eigvalsh
    # getForceMatrix() returns a numpy array... OK
    eig = LA.eigvalsh( getForceMatrix() ) # eig is a numpy array

    eig = getRidZeros(eig) # paper: "zeros do not contribute to DOS"
    # eig is now an array, not a numpy array!

    D = len(eig) # D = 3N-3, but here it's not necessary to substract -3 since
                  # we already got rid of zeros.
    Dm = D / 2.0
    N  = (D + 3) / 3.0 # comes from solving D = 3N-3
    pi = math.pi

    c  = coef(N, V)
    m1 = productSqE(eig) # not zeros!
    m2 = ( 2 * dE ) ** ( Dm  - 1 )
    m3 = 2 * ( pi ** Dm ) / math.gamma(Dm)

    return c * m1 * m2 * m3
#


# function to calculate the harmonic DOS
def harmonicDOSParticles(V, dE):
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
    rotated by an angle theta that tan(2theta) = B / (A-C)
    A particular case: V(x1,x2) = A*x1^2 - 2*sqrt(A*C) + C*x2^2
                                = (A*x1 - B*x2)^2
    which can represent the interaction of two particles joined by a string.

    """
    import math
    from scipy import linalg as LA

    # For a complex Hermitian or real symmetric matrix: eigvalsh
    # getForceMatrix() returns a numpy array... OK
    eig = LA.eigvalsh( getForceMatrix() ) # eig is a numpy array



    eig = getRidZeros(eig) # paper: "zeros do not contribute to DOS"
    # eig is now an array, not a numpy array!

#    D = len(eig) # D = 3N-3, but here it's not necessary to substract -3 since
                  # we already got rid of zeros.

    D = len(eig) - 1 # D = 3N-3, when all are particles, no external potential. Discounting one particle. All other particles will be described respect to it.

    Dm = D / 2.0
    N  = (D + 3) / 3.0 # comes from solving D = 3N-3
    pi = math.pi

    c  = coef(N, V)
    m1 = productSqE(eig) # not zeros!
    m2 = ( 2 * dE ) ** ( Dm  - 1 )
    m3 = 2 * ( pi ** Dm ) / math.gamma(Dm)

    return c * m1 * m2 * m3
#
