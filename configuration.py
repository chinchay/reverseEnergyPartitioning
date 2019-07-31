import vectormath as vmath
from vectormath import Vector3 as Vector3

# def getOctahedral():
#     a1 = Vector3(5.35391903, 7.34897441, 9.20345327)
#     a2 = Vector3(4.83487693, 9.69665956, 7.78048212)
#     a3 = Vector3(3.48505708, 8.10763956, 11.13648625)
#     a4 = Vector3(2.73813043, 7.97569725, 8.44800613)
#     a5 = Vector3(5.58070460, 9.82836905, 10.46934040)
#     a6 = Vector3(2.96585081, 10.45497144, 9.71365159)
#     listOfPositions = [a1, a2, a3, a4, a5, a6]
#     return listOfPositions

def getLJenergy(listOfPositions, sigma, epsilon):
    n = len(listOfPositions)
    e = 0
    for i in range(n):
        for j in range(n):
            if (i != j):
                e += VLJ3( listOfPositions[i], listOfPositions[j], sigma, epsilon  )
    #
    return e / 2.0 # avoid repetition
#

# for +-dot,cross vector operations:
# https://github.com/seequent/vectormath
# https://github.com/seequent/vectormath/blob/master/tests/test_vector3.py
def VLJ3(rA, rB, sigma, epsilon):
    import math
    # from vectormath import Vector3 as Vector3

    """
    Lennard-Jones potential
    eig  is the depth of the potential well.
    At rm, the potential function has the value -eig.
    """

    r = math.sqrt( (rA[0]-rB[0])**2 + (rA[1]-rB[1])**2 + (rA[2]-rB[2])**2)

    # test if atomA and atomB are in the same position!
    # delta = 1.0E-5
    # if ( (-1.0E-5 < r) and (r < 1.0E-5) ):
    #     print(rA)
    #     print(rB)
    #
    #     raise Exception('r == 0 ?? The value of r was: {}'.format(r))
    # ####

    l = sigma / r
    return 4 * epsilon * ( (l ** 12) - (l ** 6) )
#

def getLJeigenvalues(listOfPositions, epsilon, sigma,rc, getForceMatrix):
    from ase import Atoms
    from ase.build import bulk
    from ase.calculators.lj import LennardJones
    from ase.phonons import Phonons
    import numpy as np
    from scipy import linalg as LA
    # from gpaw import GPAW, FermiDirac
    # calc = LennardJones() #a.set_calculator(calc)

    # atoms = bulk('Si', 'diamond', a=5.4)
    # atoms = bulk('H', 'fcc', a=1.1, cubic=True)
    #atoms = Atoms('N3', [(0, 0, 0), (0, 0, 1.1), (0, 0, 2.2)], calculator=LennardJones() )
    # atoms = Atoms('H2', [(0, 0, 0), (0, 0, 1.12246)], calculator=LennardJones() )
    # calc = GPAW(kpts=(5, 5, 5), h=0.2, occupations=FermiDirac(0.))

    chemStr = 'H' + str(len(listOfPositions))
    calc = LennardJones(sigma=sigma, epsilon=epsilon, rc=rc)
    atoms = Atoms(chemStr, listOfPositions, calculator=calc )

    energy = atoms.get_potential_energy()

    eig = []
    if getForceMatrix:
        ph = Phonons(atoms, calc)
        ph.run()
        ph.read(acoustic=True)
        ph.clean()

        f = ph.get_force_constant()
        # f
        # f.size
        (l,m,n) = f.shape
        if l == 1:
            ff = np.reshape(f, (m,n))
        else:
             print("error")
        #
        # ff
        eig = LA.eigvalsh( ff ) # eig is a numpy array
    #
    return energy, [float("{0:.5f}".format(eig[i])) for i in range(len(eig))]
#

def getChemStr(S):
    chemStr = ''
    for i in range(len(S)):
        chemStr += S[i]
    #
    return chemStr
#

def getLJeigenvaluesB(X, S, epsilon, sigma,rc, getForceMatrix):
    from ase import Atoms
    from ase.build import bulk
    from ase.calculators.lj import LennardJones
    from ase.phonons import Phonons
    import numpy as np
    from scipy import linalg as LA
    calc = LennardJones(sigma=sigma, epsilon=epsilon, rc=rc)
    # chemStr = 'H' + str(len(X))
    # atoms = Atoms(chemStr, X, calculator=calc )
    atoms = Atoms(getChemStr(S), X, calculator=calc )
    energy = atoms.get_potential_energy()
    eig = []
    if getForceMatrix:
        ph = Phonons(atoms, calc)
        ph.run()
        ph.read(acoustic=True)
        ph.clean()

        f = ph.get_force_constant()
        (l,m,n) = f.shape
        if l == 1:
            ff = np.reshape(f, (m,n))
        else:
             print("error")
        #
        eig = LA.eigvalsh( ff ) # eig is a numpy array
    #
    return energy, [float("{0:.5f}".format(eig[i])) for i in range(len(eig))]
#
def getLJeigenvalues2B(X, S, epsilon, sigma,rc, getForceMatrix, aCell):
    from ase import Atoms
    from ase.build import bulk
    from ase.phonons import Phonons
    import numpy as np
    from scipy import linalg as LA
    from ase import Atom, Atoms
    from lammpslib import LAMMPSlib

    # chemStr = 'H' + str(len(X))
    # struct = Atoms(chemStr, X, cell=(aCell, aCell, aCell), pbc=True)
    struct = Atoms(getChemStr(S), X, cell=(aCell, aCell, aCell), pbc=True)
    lammps_header = [ "units       metal"]
    cmds          = [ "pair_style  mlip  /Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/mlip_LJ.ini",\
                      "pair_coeff  * * " ]
    mylammps = LAMMPSlib(lmpcmds = cmds, atom_types={1:1}, keep_alive=True, log_file='/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/log.txt')
    struct.set_calculator(mylammps)
    energy = struct.get_potential_energy()
    eig = []
    if getForceMatrix:
        ph = Phonons(struct, mylammps)
        ph.run()
        ph.read(acoustic=True)
        ph.clean()
        f = ph.get_force_constant()
        (l,m,n) = f.shape
        if l == 1:
            ff = np.reshape(f, (m,n))
        else:
             print("error")
        #
        eig = LA.eigvalsh( ff ) # eig is a numpy array
    #
    return energy, [float("{0:.5f}".format(eig[i])) for i in range(len(eig))]
#

def getLJeigenvalues2(listOfPositions, epsilon, sigma,rc, getForceMatrix, aCell):
    from ase import Atoms
    from ase.build import bulk
    # from ase.calculators.lj import LennardJones
    from ase.phonons import Phonons
    import numpy as np
    from scipy import linalg as LA
    # from gpaw import GPAW, FermiDirac
    # calc = LennardJones() #a.set_calculator(calc)

    # atoms = bulk('Si', 'diamond', a=5.4)
    # atoms = bulk('H', 'fcc', a=1.1, cubic=True)
    #atoms = Atoms('N3', [(0, 0, 0), (0, 0, 1.1), (0, 0, 2.2)], calculator=LennardJones() )
    # atoms = Atoms('H2', [(0, 0, 0), (0, 0, 1.12246)], calculator=LennardJones() )
    # calc = GPAW(kpts=(5, 5, 5), h=0.2, occupations=FermiDirac(0.))


    # d =  1.122 # = 2**(1/6)
    # a = 10.00
    # struct = Atoms( 'H2', positions=[(0, 0, 0), (0, 0, d)] , cell=(a, a, a), pbc=True )

    chemStr = 'H' + str(len(listOfPositions))
    # struct = Atoms(chemStr, listOfPositions, cell=(aCell, aCell, aCell)) # <<< without pbc=True you would need a very large aCell value!
    struct = Atoms(chemStr, listOfPositions, cell=(aCell, aCell, aCell), pbc=True)
    # struct = Atoms(chemStr, positions=positions , cell=(aCell, aCell, aCell), pbc=True )

    ############################################################################
    # from ase.calculators.lj import LennardJones
    # calc = LennardJones(sigma=sigma, epsilon=epsilon, rc=rc)
    # struct = Atoms(chemStr, listOfPositions, calculator=calc )

    ############################################################################
    from ase import Atom, Atoms
    from lammpslib import LAMMPSlib
    # lammps_header=['units       metal'    ,\
    #                'boundary    p p p '   ,\
    #                "atom_style	atomic"   ,\
    #                "atom_modify	map hash"    ]
    lammps_header = [ "units       metal"]
    cmds          = [ "pair_style  mlip  /Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/mlip_LJ.ini",\
                      "pair_coeff  * * " ]
    # cmds = ["pair_style    mlip  /Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/mlip_LJ.ini",\
    #         "pair_coeff    * * "       ,\
    #         "neighbor      1.5 bin "       ]

    # cmds = ["pair_style    mlip  /Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/mlip_test.ini",\
    #         "pair_coeff    * * "       ,\
    #         "neighbor      1.5 bin "       ]

    mylammps = LAMMPSlib(lmpcmds = cmds, atom_types={1:1}, keep_alive=True, log_file='/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/log.txt')
    # struct = Atoms(chemStr, listOfPositions, calculator=mylammps )
    struct.set_calculator(mylammps)
    ############################################################################


    energy = struct.get_potential_energy()

    eig = []
    if getForceMatrix:
        # ph = Phonons(struct, calc)
        ph = Phonons(struct, mylammps)
        ph.run()
        ph.read(acoustic=True)
        ph.clean()

        f = ph.get_force_constant()
        # f
        # f.size
        (l,m,n) = f.shape
        if l == 1:
            ff = np.reshape(f, (m,n))
        else:
             print("error")
        #
        # ff
        eig = LA.eigvalsh( ff ) # eig is a numpy array
    #
    return energy, [float("{0:.5f}".format(eig[i])) for i in range(len(eig))]
#


# listOfPositions = getOctahedral()
# epsilon = 0.1
# sigma   = 2.5
# rc      = 7.50
# getForceMatrix = True
# Emin, eig = getLJeigenvalues(listOfPositions, epsilon, sigma,rc, getForceMatrix)

# Emin = getLJenergy(listOfPositions)

# Emin
# eig

###############################################################################
def getRidZeros(eig):
    # harmonic oscillator has non-negative eigenvalues!
    # eig[k] != 0.0
    precision = 0.0001  # product of eigenValues will explode if we don't get rid of those "zeros" that are not zeros but ~1.0E-14
    positives = []
    [positives.append(eig[k]) if eig[k] > precision else None for k in range(len(eig)) ]
    # [positives.append(abs( eig[k] )) if abs(eig[k]) > precision else None for k in range(len(eig)) ]
    return positives
#

def productSqE(eig):
    import math
    from functools import reduce
    # http://book.pythontips.com/en/latest/map_filter.html
    product = reduce( (lambda x, y: x * y), eig )
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
    return math.factorial(m) / ( float(V) ** m ) #`float()` to ensure the division will be a float (python2 resembles Fortran in divisions)
#

# function to calculate the harmonic DOS
def harmonicDOS(dE, eig, N, V): # N=natoms, V=volume
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
    # eig = LA.eigvalsh( getForceMatrix() ) # eig is a numpy array
    # eig = LA.eigvalsh( forceMatrix ) # eig is a numpy array

    NatomsTotal = len(eig) / 3
    eig         = getRidZeros(eig) # paper: "zeros do not contribute to DOS"
    D           = len(eig) # 3N -3, or perhaps 3N - 6?
    Neffective  = D / 3
    assert(Neffective == 4)
    Dm          = D / 2.0
    c           = coef(Neffective, V)
    m1          = productSqE(eig) # not zeros!
    m2          = ( 2 * dE ) ** ( Dm  - 1 )
    # m2          = ( 2 * dE ) ** ( Dm - 4 )
    m3          = 2 * ( math.pi ** Dm ) / math.gamma(Dm)
    return c * m1 * m2 * m3
    #
    # # N = len(eig) / 3 # number of atoms
    # # Dtemp = (N - 1) * 3
    #
    # NatomsTotal = len(eig)
    # eig = getRidZeros(eig) # paper: "zeros do not contribute to DOS"
    # # eig is now an array, not a numpy array!
    #
    # D = len(eig) # D = 3N-3, but here it's not necessary to substract -3 since
    #               # we already got rid of zeros.
    #
    # # if (Dtemp != D):
    # #     print("error: Dtemp = ", Dtemp, " D = ", D)
    # #     print()
    # #
    # # assert(Dtemp == D)
    #
    #
    # Dm = D / 2.0
    # # N  = (D + 3) / 3.0 # comes from solving D = 3N-3
    #
    # pi = math.pi
    #
    # # n, nAtoms, nNeighbors, eNull, rm, eUnit, alat, A, rm = getConstants()
    # # N = nAtoms
    # # V = getVolume(A[0], A[2], A[4])
    # c = coef(N, V)
    #
    # m1 = productSqE(eig) # not zeros!
    # m2 = ( 2 * dE ) ** ( Dm  - 1 )
    # m3 = 2 * ( pi ** Dm ) / math.gamma(Dm)
    #
    # # print(eig)
    # # print("                       d = ", Dm  - 1, "f*d = ", m1 * (2**( Dm  - 1 )) * m3 )
    #
    # return c * m1 * m2 * m3
#

def Integral_harmonicDOS(dE, eig, N, V): # N=natoms, V=volume
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
    NatomsTotal = len(eig) / 3
    eig         = getRidZeros(eig) # paper: "zeros do not contribute to DOS"
    D           = len(eig) # 3N -3, or perhaps 3N - 6?
    Neffective  = D / 3
    assert(Neffective == 4)
    Dm          = D / 2.0
    c           = coef(Neffective, V)
    m1          = productSqE(eig) # not zeros!
    m2          = (2 ** ( Dm  - 1 )) * (dE ** Dm)
    # m2          = ( 2 * dE ) ** ( Dm  - 1 )
    m3          = 2 * ( math.pi ** Dm ) / math.gamma(Dm)
    return c * m1 * m2 * m3
#


def getEm(alpha_m, E_m, E_mMinus1):
    minimo = min([ 1, alpha_m ])
    rLog   = log( 1 +  (1 / minimo) ) / log(2) #<<<<<<<<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<<<<<<<
    return E_m + (  (E_m - E_mMinus1) * rLog  )
#

################################################################################
##################################################################################
from ase import Atoms
from ase.build import bulk
from ase.phonons import Phonons
import numpy as np
from scipy import linalg as LA
from ase import Atom, Atoms

class Model:
    def __init__(self, path, pot, aseStruct):
        self.path = path
        self.pot  = pot
        #
        # path should not end in '/'
        # make these variables global:
        epsilon = 0.1
        sigma   = 2.5
        rcut    = 7.50
        aCell   = 15.0
        # defining dictionary of Z:LAMMPStype
        atom_types = {}
        listZ = list(set(aseStruct.get_atomic_numbers()))
        # [ (j, i+1) for i,j in enumerate( list(set(aseStruct.get_atomic_numbers())) )]
        for i in range(len(listZ)):
            atom_types[listZ[i]] = i + 1
        #
        # path = '/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning'
        log_file      = path + "/log.txt"

        if pot == "aseLJ":
            from ase.calculators.lj import LennardJones
            self.calc = LennardJones(sigma=sigma, epsilon=epsilon, rc=rc)
        else:
            from lammpslib import LAMMPSlib
            cmds = [ "pair_style  mlip  " + path + "/" + pot, "pair_coeff  * * " ]
            # lammps_header = [ "units       metal"]  <<< it does not work in definition of mylammps. Verify code
            # print(atom_types)
            # self.mylammps
            self.calc = LAMMPSlib(lmpcmds=cmds, atom_types=atom_types, keep_alive=True, log_file=log_file)
        #
    #
    def getEnergy(self, aseStruct):
        aseStruct.set_calculator(self.calc)
        return aseStruct.get_potential_energy()
    #

    def getEnergyAndEigen(self, aseStruct):
        aseStruct.set_calculator(self.calc)
        energy = aseStruct.get_potential_energy()
        eig = []
        ph = Phonons(aseStruct, self.calc)
        ph.run()
        ph.read(acoustic=True)
        ph.clean()
        f = ph.get_force_constant()
        (l,m,n) = f.shape
        if l == 1:
            ff = np.reshape(f, (m,n))
        else:
             print("error")
        #
        eig = LA.eigvalsh( ff ) # eig is a numpy array
        #
        return energy, [float("{0:.5f}".format(eig[i])) for i in range(len(eig))]
    #
################################################################################
