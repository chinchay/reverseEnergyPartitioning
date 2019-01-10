# POTENTIAL

def VLJ(atomA, atomB):
# for +-dot,cross vector operations:
# https://github.com/seequent/vectormath
# https://github.com/seequent/vectormath/blob/master/tests/test_vector3.py
    import numpy as np
    import vectormath as vmath
    from vectormath import Vector3 as Vector3

    """
    Lennard-Jones potential
    ε  is the depth of the potential well.
    At rm, the potential function has the value −ε.
    """
    rm  = 1.0
    ε  = 1.0

    # getting atoms separation:
    rVector = atomA.position - atomB.position
    r = rVector.length # module

    # test if atomA and atomB are in the same position!
    δ = 1.0E-5
    if ( (-1.0E-5 < r) and (r < 1.0E-5) ):
        raise Exception('r == 0 ?? The value of r was: {}'.format(x))

    # Pauli repulsion at short ranges due to overlapping electron orbitals
    repulsive  = (rm / r) ** 12

    # attraction at long ranges (van der Waals force, or dispersion force??)
    attractive = -2 * ( (rm / r) ** 6 )

    return ε * ( repulsive + attractive )
#

class Atom:
    def __init__(self, name, position):
        self.name = name
        self.position = position
#



def getAtomsAtEquilibriumPositions():
    import numpy as np
    import vectormath as vmath
    from vectormath import Vector3 as Vector3
    # BCC:
    A1 = Atom("Co", Vector3(1.0, 0.0, 0.0) )
    A2 = Atom("Co", Vector3(0.0, 1.0, 0.0) )
    A3 = Atom("Co", Vector3(0.0, 0.0, 1.0) )
    A4 = Atom("Co", Vector3(0.0, 0.0, 0.0) )
    A5 = Atom("Co", Vector3(1.0, 1.0, 0.0) )
    A6 = Atom("Co", Vector3(1.0, 0.0, 1.0) )
    A7 = Atom("Co", Vector3(0.0, 1.0, 1.0) )
    A8 = Atom("Co", Vector3(1.0, 1.0, 1.0) )
    #lAtoms = [A1,A2,A3,A4,A5,A6,A7,A8]
    lAtoms = [A1,A2,A3]
    return lAtoms
#


def getTotalPotentialEnergy(lAtoms, V): # V is a function name `VLJ` for example
    n = len(lAtoms)
    pot = 0.0
    for i in range(n):
        for j in range(n):
            if (j > i): # avoid repeated interactions, and self-interaction
                pot += V( lAtoms[i], lAtoms[j] )
    #
    return pot
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
