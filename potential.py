# POTENTIAL

def VLJ(r):
# for +-dot,cross vector operations:
# https://github.com/seequent/vectormath
# https://github.com/seequent/vectormath/blob/master/tests/test_vector3.py
    import numpy as np
    import vectormath as vmath
    # from vectormath import Vector3 as Vector3

    """
    Lennard-Jones potential
    eps  is the depth of the potential well.
    At rm, the potential function has the value -eps.
    """
    rm  = 1.0
    eps  = 1.0

    # getting atoms separation:
    # rVector = atomA.position - atomB.position
    # r = rVector.length # module

    # test if atomA and atomB are in the same position!
    delta = 1.0E-5
    if ( (-1.0E-5 < r) and (r < 1.0E-5) ):
        raise Exception('r == 0 ?? The value of r was: {}'.format(x))

    # Pauli repulsion at short ranges due to overlapping electron orbitals
    repulsive  = (rm / r) ** 12

    # attraction at long ranges (van der Waals force, or dispersion force??)
    attractive = -2 * ( (rm / r) ** 6 )

    return eps * ( repulsive + attractive )
#

def VLJ2(rA, rB):
# for +-dot,cross vector operations:
# https://github.com/seequent/vectormath
# https://github.com/seequent/vectormath/blob/master/tests/test_vector3.py
    import numpy as np
    import vectormath as vmath
    import math
    # from vectormath import Vector3 as Vector3

    """
    Lennard-Jones potential
    eps  is the depth of the potential well.
    At rm, the potential function has the value -eps.
    """
    rm  = 1.0
    eps  = 1.0

    # getting atoms separation:
    rVector = [rA[i] - rB[i] for i in range(3)]
    r = math.sqrt(rVector[0]**2 + rVector[1]**2 + rVector[2]**2)

    # test if atomA and atomB are in the same position!
    delta = 1.0E-5
    if ( (-1.0E-5 < r) and (r < 1.0E-5) ):
        raise Exception('r == 0 ?? The value of r was: {}'.format(x))

    # Pauli repulsion at short ranges due to overlapping electron orbitals
    repulsive  = (rm / r) ** 12

    # attraction at long ranges (van der Waals force, or dispersion force??)
    attractive = -2 * ( (rm / r) ** 6 )

    # cutoff:
    # it avoids particles on the diagonal sqrt(2)*lattice, only first
    # nearest neighbors.
    # if (r > 1.3 * rm):
        # return 0

    return eps * ( repulsive + attractive )
#

class Atom:
    def __init__(self, name, position, isSlavePeriodic, indexMaster):
        self.name = name
        self.position = position
        self.isSlavePeriodic = isSlavePeriodic # if atoms is subject to copy position of the master atom at the border (periodicity)
        self.indexMaster = indexMaster # if not isSlavePeriodic -> withWhomPeriodic==-1
        # Make sure indexSlave < indexMaster  !!!
#



def getAtomsAtEquilibriumPositions():
    import numpy as np
    from math import sqrt, cos, sin, pi
    import vectormath as vmath
    from vectormath import Vector3 as Vector3
    # BCC:
    # A1 = Atom("Co", Vector3(1.0, 0.0, 0.0) )
    # A2 = Atom("Co", Vector3(0.0, 1.0, 0.0) )
    # A3 = Atom("Co", Vector3(0.0, 0.0, 1.0) )
    # A4 = Atom("Co", Vector3(0.0, 0.0, 0.0) )
    # A5 = Atom("Co", Vector3(1.0, 1.0, 0.0) )
    # A6 = Atom("Co", Vector3(1.0, 0.0, 1.0) )
    # A7 = Atom("Co", Vector3(0.0, 1.0, 1.0) )
    # A8 = Atom("Co", Vector3(1.0, 1.0, 1.0) )





    a = sqrt(2.0)
    n = 3 # for a cuboid of n^3 atoms
    ex = Vector3(1, 0, 0)
    ey = Vector3(0, 1, 0)
    ez = Vector3(0, 0, 1)
    a1 = (a / 2) * (ex + ey)
    a2 = (a / 2) * (ey + ez)
    a3 = (a / 2) * (ex + ez)

    lAtoms = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                vx = i * a1
                vy = j * a2
                vz = k * a3
                lAtoms.append( Atom("H", vx + vy + vz, False, -1 ) )





    # a = 0.978351 #0.982
    # nMax = 5 #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #
    # lAtoms = []
    # for i in range(nMax):
    #     for j in range(nMax):
    #         # k = (i * nMax) + j
    #         x = i * a
    #         y = j * a
    #         lAtoms.append( Atom("H", Vector3(x, y, 0.0), False, -1 ) )
    #



    #
    # A0 = Atom("Co", Vector3(0, 0, 0), False, -1 )
    #
    # # dist = a
    # A1 = Atom("Co", Vector3(a, 0,  0), False, -1 )
    # A2 = Atom("Co", Vector3(-a, 0, 0), False, -1 )
    # A3 = Atom("Co", Vector3(0, a,  0), False, -1 )
    # A4 = Atom("Co", Vector3(0, -a, 0), False, -1 )
    #
    # # dist = a*sqrt(2)
    # lx = a
    # ly = a
    # A5 = Atom("Co", Vector3(lx,  ly, 0), False, -1 )
    # A6 = Atom("Co", Vector3(lx, -ly, 0), False, -1 )
    # A7 = Atom("Co", Vector3(-lx, ly, 0), False, -1 )
    # A8 = Atom("Co", Vector3(-lx, -ly,0), False, -1 )
    #
    # # dist = 2a
    # l = 2 * a
    # A9  = Atom("Co", Vector3( l, 0, 0), False, -1 )  #1
    # A10 = Atom("Co", Vector3(-l, 0, 0), True,   9 )  #1
    # A11 = Atom("Co", Vector3(0,  l, 0), False, -1 )  #2
    # A12 = Atom("Co", Vector3(0, -l, 0), True,  11 )  #2
    #
    # # dist = a*sqrt(5)
    # lx = a
    # ly = 2 * a
    # A13 = Atom("Co", Vector3(lx,  ly, 0), False, -1 )  #3
    # A14 = Atom("Co", Vector3(lx, -ly, 0), True,  13 )  #3
    # A15 = Atom("Co", Vector3(-lx, ly, 0), False, -1 )  #4
    # A16 = Atom("Co", Vector3(-lx, -ly,0), True,  13 )  #4
    #
    # lx = 2 * a
    # ly = a
    # A17 = Atom("Co", Vector3(lx,  ly, 0), False, -1 )  #5
    # A18 = Atom("Co", Vector3(lx, -ly, 0), False, -1 )  #6
    # A19 = Atom("Co", Vector3(-lx, ly, 0), True,  18 )  #5
    # A20 = Atom("Co", Vector3(-lx, -ly,0), True,  17 )  #6
    #
    # lx = 2 * a
    # ly = 2 * a
    # A21 = Atom("Co", Vector3(lx,  ly, 0), False, -1 )  #7
    # A22 = Atom("Co", Vector3(lx, -ly, 0), True,  21 )  #7
    # A23 = Atom("Co", Vector3(-lx, ly, 0), True,  21 )  #7
    # A24 = Atom("Co", Vector3(-lx, -ly,0), True,  21 )  #7
    #
    #
    # lAtoms = [A0, A1,A2,A3,A4,  A5,A6,A7,A8,A9,A10,\
    #           A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,\
    #           A21,A22,A23,A24]
    #
    #
    #
    # # lAtoms = [A1,A2,A3,A4,A5,A6,A7,A8]
    # # lAtoms = [A1,A2,A3]
    return lAtoms, a1, a2, a3
#

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

def get3DPotNearFirstNeighb(lAtoms, V, a1, a2, a3):
    from vectormath import Vector3 as Vector3
    # from math import sqrt
    n = int(  (float(len(lAtoms)) +0.1) ** (1./3)  ) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # print(len(lAtoms), n)

    A1 = n * a1
    A2 = n * a2
    A3 = n * a3

    pot = 0
    for i in range(n):
        # find index for right and left neighbors using periodic boundary conditions
        iplus1, imins1, iRp1, iRm1 = getNeighborsIndx(i, n)
        for j in range(n):
            # find index for top and bottom neighbors using periodic boundary conditions
            jplus1, jmins1, jRp1, jRm1 = getNeighborsIndx(j, n)
            for k in range(n):
                # find index for up and down neighbors using periodic boundary conditions
                kplus1, kmins1, kRp1, kRm1 = getNeighborsIndx(k, n)

                # Energy between nearest first neighbors:
                lc = (i * n * n) + (j * n) + k
                l1 = (iplus1 * n * n) + (j * n) + k
                l2 = (imins1 * n * n) + (j * n) + k
                l3 = (i * n * n) + (jplus1 * n) + k
                l4 = (i * n * n) + (jmins1 * n) + k
                l5 = (i * n * n) + (j * n) + kplus1
                l6 = (i * n * n) + (j * n) + kmins1

                rc = lAtoms[lc].position
                r1 = lAtoms[l1].position + (iRp1 * A1)
                r2 = lAtoms[l2].position + (iRm1 * A1)
                r3 = lAtoms[l3].position + (jRp1 * A2)
                r4 = lAtoms[l4].position + (jRm1 * A2)
                r5 = lAtoms[l5].position + (kRp1 * A3)
                r6 = lAtoms[l6].position + (kRm1 * A3)

                pot += V( (rc - r1).length  ) +\
                       V( (rc - r2).length  ) +\
                       V( (rc - r3).length  ) +\
                       V( (rc - r4).length  ) +\
                       V( (rc - r5).length  ) +\
                       V( (rc - r6).length  )
    #
    return pot / 2  # Because of double counting.

####
def getTotalPotentialEnergy(lAtoms, V): # V is a function name `VLJ` for example
    n = len(lAtoms)
    pot = 0.0
    i = 0

    from math import sqrt
    nMax = int(sqrt(n)) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # for each step in time, loop over every cell and determine flip

    pot = 0
    for i in range(nMax):

        # find index for right and left neighbors using periodic boundary conditions
        iplus1 = i + 1
        if i == nMax-1: iplus1 = 0
        imins1 = i - 1
        if i == 0: imins1 = nMax-1

        for j in range(nMax):

            # find index for top and bottom neighbors using periodic boundary conditions
            jplus1 = j + 1
            if j == nMax-1: jplus1 = 0
            jmins1 = j - 1
            if j == 0: jmins1 = nMax-1

            # Energy between nearest first neighbors:
            kc  = (i * nMax) + j
            k1  = (iplus1 * nMax) + j
            k2  = (imins1 * nMax) + j
            k3  = (i * nMax) + jplus1
            k4  = (i * nMax) + jmins1
            pot += V( lAtoms[kc], lAtoms[k1]  ) +\
                   V( lAtoms[kc], lAtoms[k2]  ) +\
                   V( lAtoms[kc], lAtoms[k3]  ) +\
                   V( lAtoms[kc], lAtoms[k4]  )
        #
    #

    #
    # for i in range(n):
    #     for j in range(n):
    #         if (j != i): # avoid repeated interactions, and self-interaction
    #             pot += V( lAtoms[i], lAtoms[j] )
    #
    return pot / 2 # Because of double counting.
#


def test():
# for +-dot,cross vector operations:
# https://github.com/seequent/vectormath
# https://github.com/seequent/vectormath/blob/master/tests/test_vector3.py
    import numpy as np
    import vectormath as vmath
    from vectormath import Vector3 as Vector3
    A = Atom("Co", Vector3(0.0, 1.0, 0.0) )
    B = Atom("Co", Vector3(1.0, 0.0, 0.0) )
    C = Atom("Co", Vector3(0.0, 1.0, 3.0) )
    assert VLJ(A, B) == -0.2343749999999999
    assert VLJ(B, C) == -0.001502065127873102


#


# # for +-dot,cross vector operations:
# # https://github.com/seequent/vectormath
# # https://github.com/seequent/vectormath/blob/master/tests/test_vector3.py
# import numpy as np
# import vectormath as vmath
# from vectormath import Vector3 as Vector3
#
# v = Vector3(1,0,0)
# a = Vector3(0,1,0)
# v
# v * a
# v.cross(a)
# v.dot(v)
# a.normalize()
# a = Vector3(3,4,0)
# a.length
# a = [1,2,3,0,9]
# len(a)
