################################################################################
# python3
import os
print(os.getcwd())
os.chdir('/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/')
print(os.getcwd())

# import graphlab
################################################################################

import vectormath as vmath
from vectormath import Vector3 as Vector3


def kroneckerDelta(i, j):
    if (i == j):
        return 1
    return 0
#

def getVector(i, j, k, alat):
    return (i * alat[0]) + (j * alat[1]) + (k * alat[2])
#

def getIndexAtom(i, j, k, n):
    return (i * n * n) + (j * n) + k

def dist(r, rCenter, dR):
    d = r + dR - rCenter
    return d.length

def getNeighborsIndx(i, n):
    iplus1 = i + 1
    iRp1 = 0
    if i == n - 1:
        iplus1 = 0
        iRp1 = 1 # translate the neighbor's position to consider periodicity
    #
    imins1 = i - 1
    iRm1 = 0
    if i == 0:
        imins1 = n - 1
        iRm1 = -1
    #
    return iplus1, imins1, iRp1, iRm1

#

def getConstants():
    from math import sqrt
    n = 3 #3
    nAtoms = n ** 3
    nNeighbors = 6
    rm = 1

    eNull = Vector3(0, 0, 0)
    rm = 1
    ex = Vector3(1, 0, 0)
    ey = Vector3(0, 1, 0)
    ez = Vector3(0, 0, 1)
    sq2s2 = sqrt(2.0)/2

    eUnit = [ex, ey, ez]
    alat =  [sq2s2 * (ex + ey), sq2s2 * (ey + ez), sq2s2 * (ex + ez)]

    A = [ n * alat[0],\
          n * alat[0],\
          n * alat[1],\
          n * alat[1],\
          n * alat[2],\
          n * alat[2]  ]

    return n, nAtoms, nNeighbors, eNull, rm, eUnit, alat, A, rm


def addIfNotCounted(iNeighbor, rc, r, lc, l, dR, countedInPot):
    from potential import VLJ as V
    dPot = 0
    if ( (countedInPot[lc][l] == False) or (countedInPot[l][lc]) == False):
        dPot = V( (r + dR[lc][iNeighbor-1] - rc).length  )
        countedInPot[lc][l] == True
        countedInPot[l][lc] == True
    return dPot

def loop(positionAtomsInput, isFirstTime):
    import potential
    from potential import VLJ as V

    n, nAtoms, nNeighbors, eNull, rm, eUnit, alat, A, rm = getConstants()

    positionAtoms = []
    posNeigh = []
    neb = []
    energyPerParticle = []

    countedInPot = [ [ False for i in range(nAtoms) ] for j in range(nAtoms) ]

    areNeb = [ [ False for i in range(nAtoms) ] for j in range(nAtoms) ]
    # areNeb

    d  = [0 for i in range(nNeighbors)]
    dR = [ [ eNull for i in range(nNeighbors) ] for j in range(nAtoms) ]
    # dR

    pot = 0
    contador = 0
    for i in range(n):
        # find index for right and left neighbors using periodic boundary conditions
        iplus1, imins1, iRp1, iRm1 = getNeighborsIndx(i, n)
        for j in range(n):
            # find index for top and bottom neighbors using periodic boundary conditions
            jplus1, jmins1, jRp1, jRm1 = getNeighborsIndx(j, n)
            for k in range(n):
                # find index for up and down neighbors using periodic boundary conditions
                kplus1, kmins1, kRp1, kRm1 = getNeighborsIndx(k, n)
                contador += 1

                # Energy between nearest first neighbors:
                l0 = getIndexAtom(i,      j,      k,      n)
                l1 = getIndexAtom(iplus1, j,      k,      n)
                l2 = getIndexAtom(imins1, j,      k,      n)
                l3 = getIndexAtom(i,      jplus1, k,      n)
                l4 = getIndexAtom(i,      jmins1, k,      n)
                l5 = getIndexAtom(i,      j,      kplus1, n)
                l6 = getIndexAtom(i,      j,      kmins1, n)

                if isFirstTime: # positions are not defined, return positions at their equilibrium state:
                    r0 = getVector(i,      j,      k,      alat)
                    r1 = getVector(iplus1, j,      k,      alat)
                    r2 = getVector(imins1, j,      k,      alat)
                    r3 = getVector(i,      jplus1, k,      alat)
                    r4 = getVector(i,      jmins1, k,      alat)
                    r5 = getVector(i,      j,      kplus1, alat)
                    r6 = getVector(i,      j,      kmins1, alat)

                    out = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/equilibriumPositions.txt"
                    if contador == 1:
                        f = open(out, "w")
                    else:
                        f = open(out, "a")
                    #
                    f.write(  str(r0[0]) + "   " + str(r0[1]) + "   " + str(r0[2]) + "\n" )
                    f.close()

                else: # update positions to calculate energy of lattice vibrating:
                    r0 = positionAtomsInput[l0]
                    r1 = positionAtomsInput[l1]
                    r2 = positionAtomsInput[l2]
                    r3 = positionAtomsInput[l3]
                    r4 = positionAtomsInput[l4]
                    r5 = positionAtomsInput[l5]
                    r6 = positionAtomsInput[l6]
                ##

                positionAtoms.append(r0)

                posNeigh.append( [l0, r1, r2, r3, r4, r5, r6] )
                neb.append( [l0, l1, l2, l3, l4, l5, l6] )

                # l1
                # l0
                # areNeb[l1][l0]

                areNeb[l1][l0] = True
                areNeb[l2][l0] = True
                areNeb[l3][l0] = True
                areNeb[l4][l0] = True
                areNeb[l5][l0] = True
                areNeb[l6][l0] = True

                d[0] =  iRp1
                d[1] =  iRm1
                d[2] =  jRp1
                d[3] =  jRm1
                d[4] =  kRp1
                d[5] =  kRm1

                for m in range(nNeighbors):
                    dR[l0][m] = d[m] * A[m]
                #


                # energyParticle = V( (r1 + dR[l0][0] - r0).length  ) +\
                #                  V( (r2 + dR[l0][1] - r0).length  ) +\
                #                  V( (r3 + dR[l0][2] - r0).length  ) +\
                #                  V( (r4 + dR[l0][3] - r0).length  ) +\
                #                  V( (r5 + dR[l0][4] - r0).length  ) +\
                #                  V( (r6 + dR[l0][5] - r0).length  )
                # energyParticle = energyParticle / 2
                # energyPerParticle.append(energyParticle)
                #
                # pot += energyParticle


                # addIfNotCounted(iNeighbor, rc, r, lc, l, dR, countedInPot)

                pot += addIfNotCounted(1, r0, r1, l0, l1, dR, countedInPot)
                pot += addIfNotCounted(2, r0, r2, l0, l2, dR, countedInPot)
                pot += addIfNotCounted(3, r0, r3, l0, l3, dR, countedInPot)
                pot += addIfNotCounted(4, r0, r4, l0, l4, dR, countedInPot)
                pot += addIfNotCounted(5, r0, r5, l0, l5, dR, countedInPot)
                pot += addIfNotCounted(6, r0, r6, l0, l6, dR, countedInPot)


                # pot += V( (r1 + dR[l0][0] - r0).length  ) +\
                #        V( (r2 + dR[l0][1] - r0).length  ) +\
                #        V( (r3 + dR[l0][2] - r0).length  ) +\
                #        V( (r4 + dR[l0][3] - r0).length  ) +\
                #        V( (r5 + dR[l0][4] - r0).length  ) +\
                #        V( (r6 + dR[l0][5] - r0).length  )

    #
    # Divide by 6:
    # pot = pot / 6
    # pot = pot / 2
    return pot, energyPerParticle, positionAtoms, posNeigh, neb, areNeb, dR, n, nAtoms, nNeighbors, rm, areNeb, neb    # Because of double counting.

####

def cosine(l, p, i, r0, posNeigh, dR, x):
    # L = l + 1 # index neighbors begins in zero
    # cosine = (posNeigh[p][l+1][i] + dR[p][l][i] - x[p][i]) / r0
    # return (cosine / 6) / 2  # <<<<< 6 is the number of neighbors in 3D

    # absdR = abs(dR[p][l][i])
    # if ( ( 0.8 < absdR) and (absdR < 1.2) ):
    #     return 0
    # else:
    #     return (posNeigh[p][l+1][i] + dR[p][l][i] - x[p][i]) / r0
    #


    return (posNeigh[p][l+1][i] + dR[p][l][i] - x[p][i]) / r0
    ## c[l_, p_, i_, r0_] := ( (posNeigh[[p, l, i]] + dR[[p, l, i]])  -   x[[p, i]])/r0;

def h1(q, j, p, i,posNeigh, dR, x, nNeighbors, rm):
    if p == q:
        sum = 0
        for l in range(nNeighbors):
            sum += cosine(l, p, j, rm, posNeigh, dR, x) * cosine(l, p, i, rm, posNeigh, dR, x)
        #
        return sum
    #
    return 0
#

def h2(q, j, p, i,posNeigh, dR, x, nNeighbors, rm, areNeb, neb):
    if ((p != q) and areNeb[p][q]):
        sum = 0
        for l in range(nNeighbors):
            L = l + 1 # index neighbors begins in zero
            sum += -cosine(l, p, j, rm,posNeigh, dR, x) * cosine(l, p, i, rm,posNeigh, dR, x) * kroneckerDelta(neb[p][L], q)
        #
        return sum
    #
    return 0
#

def hSum(q, j, p, i,posNeigh, dR, x, nNeighbors, rm, areNeb, neb):

    return h1(q, j, p, i,posNeigh, dR, x, nNeighbors, rm) + h2(q, j, p, i,posNeigh, dR, x, nNeighbors, rm, areNeb, neb)
#

def indxP(s):
    return s // 3

def indxI(s):
    return s % 3

def force(s, t,posNeigh, dR, x, nNeighbors, rm, areNeb, neb):

    return hSum(indxP(s), indxI(s), indxP(t), indxI(t),posNeigh, dR, x, nNeighbors, rm, areNeb, neb)
#


def getConfigEquilibrium():
    import copy

    positionAtomsInput = []
    isFirstTime = True
    potE, energyPerParticle, positionAtoms, posNeigh, neb, areNeb, dR, n, nAtoms, nNeighbors, rm, areNeb, neb = loop(positionAtomsInput, isFirstTime)
    n, nAtoms, nNeighbors, eNull, rm, eUnit, alat, A, rm = getConstants()
    Emin = potE

    Xeq = copy.deepcopy(positionAtoms)
    a1 = alat[0]
    a2 = alat[1]
    a3 = alat[2]

    ######## forceMatrix
    rm = 1
    nOrder = n
    # x = copy.deepcopy(positionAtoms)

    # factor 72 is due to the second radial derivative of Lennard-Jones potential: V''(r) = 72 evaluated at rm = 1 for every particle in the crystal (equilibrium)
    # factor * 2 due to correction in the calculations to find the forceMatrix

    fact = 72 #72 / 2
    fact = fact * 2
    forceMatrix = [  [fact * force(s, t, posNeigh, dR, Xeq, nNeighbors, rm, areNeb, neb) for s in range(3 * nAtoms)] for t in range(3 * nAtoms)]
    # temp = [  [fact * force(s, t, posNeigh, dR, Xeq, nNeighbors, rm, areNeb, neb) for s in range(3 * nAtoms)] for t in range(3 * nAtoms)]
    # forceMatrix = copy.deepcopy(temp)
    #################

    return Xeq, a1, a2, a3, Emin, forceMatrix

#
#
# def getForceMatrix():
#     import copy
#
#     positionAtomsInput = []
#     isFirstTime = True
#     potE, energyPerParticle, positionAtoms, posNeigh, neb, areNeb, dR, n, nAtoms, nNeighbors, rm, areNeb, neb = loop(positionAtomsInput, isFirstTime)
#
#     rm = 1
#     nOrder = n
#
#     x = copy.deepcopy(positionAtoms)
#
#     forceMatrix = [ [force(s, t, posNeigh, dR, x, nNeighbors, rm, areNeb, neb) for s in range(3 * nAtoms)] for t in range(3 * nAtoms)]
#     return forceMatrix
# ####

# def getAtomsAtEquilibriumPositions():
#     import copy
#
#     positionAtomsInput = []
#     isFirstTime = True
#     potE, energyPerParticle, positionAtoms, posNeigh, neb, areNeb, dR, n, nAtoms, nNeighbors, rm, areNeb, neb = loop(positionAtomsInput, isFirstTime)
#     n, nAtoms, nNeighbors, eNull, rm, eUnit, alat, A, rm = getConstants()
#
#     x = copy.deepcopy(positionAtoms)
#     a1 = alat[0]
#     a2 = alat[1]
#     a3 = alat[2]
#     return x, a1, a2, a3, potE
# #

################################################################################
def getRidZeros(ε):
    # harmonic oscillator has non-negative eigenvalues!
    # ε[k] != 0.0
    precision = 0.0001  # product of eigenValues will explode if we don't get rid of those "zeros" that are not zeros but ~1.0E-14
    positives = []
    [positives.append(ε[k]) if ε[k] > precision else None for k in range(len(ε)) ]
    return positives
#

def productSqE(ε):
    import math
    from functools import reduce
    # http://book.pythontips.com/en/latest/map_filter.html
    product = reduce( (lambda x, y: x * y), ε )
    return 1 / math.sqrt( product )
#

def getVolume(A, B, C):
    import numpy as np

    crossProd = np.cross(A, B)
    return np.dot(crossProd, C)

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

# def getProdInvSqrtEigenvals():
#     from numpy import linalg as LA
#
#     # For a complex Hermitian or real symmetric matrix: eigvalsh
#     # getForceMatrix() returns a numpy array... OK
#     ε = LA.eigvalsh( getForceMatrix() ) # ε is a numpy array
#
#     ε = getRidZeros(ε) # paper: "zeros do not contribute to DOS" <<< still wondering why so many zeros???????????
#     # ε is now an array, not a numpy array!
#     return productSqE(ε)

# function to calculate the harmonic DOS
def harmonicDOS(dE, forceMatrix):
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
    from scipy import linalg as LA

    # For a complex Hermitian or real symmetric matrix: eigvalsh
    # getForceMatrix() returns a numpy array... OK
    # ε = LA.eigvalsh( getForceMatrix() ) # ε is a numpy array
    ε = LA.eigvalsh( forceMatrix ) # ε is a numpy array

    ε = getRidZeros(ε) # paper: "zeros do not contribute to DOS"
    # ε is now an array, not a numpy array!

    D = len(ε) # D = 3N-3, but here it's not necessary to substract -3 since
                  # we already got rid of zeros.
    Dm = D / 2.0
    N  = (D + 3) / 3.0 # comes from solving D = 3N-3
    π = math.pi

    n, nAtoms, nNeighbors, eNull, rm, eUnit, alat, A, rm = getConstants()
    N = nAtoms
    V = getVolume(A[0], A[2], A[4])
    c  = coef(N, V)

    m1 = productSqE(ε) # not zeros!
    m2 = ( 2 * dE ) ** ( Dm  - 1 )
    m3 = 2 * ( π ** Dm ) / math.gamma(Dm)

    # print(ε)
    # print("                       d = ", Dm  - 1, "f*d = ", m1 * (2**( Dm  - 1 )) * m3 )

    return c * m1 * m2 * m3
#



################################################################################


def get3DPotNearFirstNeighb(lAtoms, V, a1, a2, a3):
    import copy

    positionAtomsInput = copy.deepcopy(lAtoms)
    isFirstTime = False
    potE, energyPerParticle, positionAtoms, posNeigh, neb, areNeb, dR, n, nAtoms, nNeighbors, rm, areNeb, neb = loop(positionAtomsInput, isFirstTime)

    return potE, energyPerParticle


def getHarmonicEnergy(X, Xeq, forceMatrix):
    import numpy as np

    eTotHarmonic = 0
    nAtoms = len(X)

    dw = []
    for i in range(nAtoms):

        # dX.append(X[i] - Xeq[i])
        dw.append(X[i][0] - Xeq[i][0])
        dw.append(X[i][1] - Xeq[i][1])
        dw.append(X[i][2] - Xeq[i][2])

        # # Do not consider translation, then describe everyone respect to one atom:
        # # dX.append(X[i] - Xeq[i])
        # dw.append( (X[i][0] - X[0][0]) - (Xeq[i][0] - Xeq[0][0]) )
        # dw.append( (X[i][1] - X[0][1]) - (Xeq[i][1] - Xeq[0][1]) )
        # dw.append( (X[i][2] - X[0][2]) - (Xeq[i][2] - Xeq[0][2]) )


    #


    deltaEharmonic = 0
    dimension = 3 * nAtoms
    for i in range(dimension):
        for j in range(dimension):
            deltaEharmonic += dw[i] * forceMatrix[i][j] * dw[j]
            # if i==j: deltaEharmonic += dw[i] * dw[j]
    #
    deltaEharmonic = deltaEharmonic / 2

    # rHyper = 0
    # for i in range(nAtoms):
    #     rHyper += X[i][0] ** 2
    #     rHyper += X[i][1] ** 2
    #     rHyper += X[i][2] ** 2
    #

    rHyper = 0
    # for i in range(dimension):
    #     rHyper += dw[i] ** 2

    for i in range(nAtoms):
        rHyper += (X[i][0]-X[0][0])**2 +\
                  (X[i][1]-X[0][1])**2 +\
                  (X[i][2]-X[0][2])**2        
    #

    rHyper = rHyper ** 0.5

    return rHyper, deltaEharmonic #+Emin
#



################################################################################
# harmonicDOS(8)
#
# n, nAtoms, nNeighbors, eNull, rm, eUnit, alat, A, rm = getConstants()
# N = nAtoms
# V = getVolume(A[0], A[2], A[4])
# coef(N, V)
#
# ε = LA.eigvalsh( getForceMatrix() )
# ε = getRidZeros(ε)
# m1 = productSqE(ε)
# m1
#
# D = len(ε)
# D
# Dm = D / 2.0
# N  = (D + 3) / 3.0 # comes from solving D = 3N-3
# π = math.pi
# m2 = ( 2 * 1 ) ** ( Dm  - 1 )
# m3 = 2 * ( π ** Dm ) / math.gamma(Dm)
#
# Dm
# m2
# 2**27
# m3
#
# f = getForceMatrix()
# f
# f[0]
#
#
# from numpy import linalg as LA
# # from math import abs
# w, v = LA.eigh( f )
#
# roundEigenvals = []
# for i in range(len(w)):
#     if  abs(w[i])> 0.0001:
#         roundEigenvals.append(w[i])
#     else:
#         roundEigenvals.append(0)
# #
# roundEigenvals
