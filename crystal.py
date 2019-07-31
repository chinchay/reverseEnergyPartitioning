################################################################################
# python3
import os
print(os.getcwd())
os.chdir('/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/')
print(os.getcwd())

# import graphlab
################################################################################

import numpy as np
import math
from pymatgen import Lattice, Structure, Molecule
from potential import VLJ as V


def getFccEquilibrium():
    # alat = 1.0
    alat = math.sqrt(2.0)
    fcc = Structure( Lattice.cubic(alat), ["Co", "Co", "Co", "Co"],
                       [[0.0, 0.0, 0.0],
                        [0.0, 0.5, 0.5],
                        [0.5, 0.0, 0.5],
                        [0.5, 0.5, 0.0]]
                      )




    # fcc = Structure( Lattice.cubic(alat), ["Co"],
    #                    [[0.0, 0.0, 0.0]]
    #                   )

    # Find the primitive unit cell:
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    sym_finder = SpacegroupAnalyzer(fcc)
    fcc = sym_finder.get_primitive_standard_structure()
    # new_structure.to(fmt='poscar', filename='POSCAR_test_2')



    # Make a supercell
    n  = 3 #3
    fcc.make_supercell([n, n, n])
    return fcc
#
#
# fcc = getFccEquilibrium()
# fcc


#
# def getNeighborsIndxs(i, struct, radio):
#     # see https://pymatgen.org/pymatgen.core.structure.html for more info
#     # determine the neighbors of each site:
#     site = struct[i]
#     tupleNeigh = struct.get_neighbors(site, radio, include_index=True)
#
#     n = len(tupleNeigh)
#     jNeighbors = [tupleNeigh[l][2] for l in range(n)]
#     diffRvec   = [tupleNeigh[l][0].coords - site.coords for l in range(n)]
#     return jNeighbors, diffRvec
# #
fcc = getFccEquilibrium()
fcc[0].distance(fcc[1])


def getNeighborsfromI(i, struct, radio):
    # see https://pymatgen.org/pymatgen.core.structure.html for more info
    # determine the neighbors of each site:
    site = struct[i]
    tupleNeigh = struct.get_neighbors(site, radio, include_index=True)

    n = len(tupleNeigh)
    indxNeigsForI = [tupleNeigh[l][2] for l in range(n)]
    dRneighborsfromI = [tupleNeigh[l][0].coords - site.coords for l in range(n)]
    return indxNeigsForI, dRneighborsfromI
#
def getDRneighborsFromEachSite(struct, radio):
    # indxNeigsFromEachSite   = []
    dRneighborsFromEachSite = []

    for i in range(struct.num_sites):
        indxNeigsForI, dRneighborsfromI = getNeighborsfromI(i, struct, radio)
        dRneighborsFromEachSite.append(dRneighborsfromI)
        # indxNeigsFromEachSite.append(indxNeigsForI)
    #
    return dRneighborsFromEachSite
#
def getIndxNeigsFromEachSite(struct, radio):
    indxNeigsFromEachSite   = []
    # dRneighborsFromEachSite = []

    for i in range(struct.num_sites):
        indxNeigsForI, dRneighborsfromI = getNeighborsfromI(i, struct, radio)
        # dRneighborsFromEachSite.append(dRneighborsfromI)
        indxNeigsFromEachSite.append(indxNeigsForI)
    #
    return indxNeigsFromEachSite
#
def getEnergy(struct, radio, dRneighborsFromEachSite): #radio=1.5
    energy = 0

    # use this when atoms have departed and others have came:
    if (len(dRneighborsFromEachSite) == 0):
        dRneighborsFromEachSite = getDRneighborsFromEachSite(struct, radio)
    #

    for i in range(struct.num_sites):
        for j in range(len(dRneighborsFromEachSite[i])):
            dR = dRneighborsFromEachSite[i][j]
            energy += V( np.linalg.norm(dR) )
    #
    return energy / 2.0  #struct.num_sites
#
def getIndxSiteOfNeighbor(m, lth_neighborOfm, struct, radio):
    indxNeigsForI, dRneighborsfromI = getDRneighborsfromI(i, struct, radio)
#
def getCos(m, lth_neighborOfm, p, i, indxNeigsFromEachSite, dRneighborsFromEachSite):
    # (X_n_i - X_m_i) / Rnm
    indxNeigSite = indxNeigsFromEachSite[m][lth_neighborOfm]
    if (p == m):
        Rvec = dRneighborsFromEachSite[m][lth_neighborOfm]
        return Rvec[i] / np.linalg.norm(Rvec)
    elif (p == indxNeigSite):
        Rvec = dRneighborsFromEachSite[m][lth_neighborOfm]
        return -Rvec[i] / np.linalg.norm(Rvec)

    #
    return 0
#
#
# len(fcc.get_neighbors(fcc[26], 1.1, include_index=True))
# indxNeigsFromEachSite[0]
# indxNeigsFromEachSite[26]
#
# i

# getCos(0, 8, 0, 1, getDRneighborsFromEachSite(fcc, 1.5))

#
# def getSumCos(m, p, i, dRneighborsFromEachSite):
#     sumCos = 0
#     for ith_neighborOfm in range(len(dRneighborsFromEachSite)):
#         sumCos += getCos(m, ith_neighborOfm, p, i, dRneighborsFromEachSite)
#     #
#     return sumCos
# #
def indxP(s):
    return s // 3
#
def indxI(s):
    return s % 3
#
def forceMatrix_st(s, t, num_sites, indxNeigsFromEachSite, dRneighborsFromEachSite):
    h = 0
    p, i, q, j = indxP(s), indxI(s), indxP(t), indxI(t)

    for m in range(num_sites):
        for lth_neighborOfm in range(len(dRneighborsFromEachSite[m])):
            cos_mli = getCos(m, lth_neighborOfm, p, i, indxNeigsFromEachSite, dRneighborsFromEachSite)
            cos_mlj = getCos(m, lth_neighborOfm, q, j, indxNeigsFromEachSite, dRneighborsFromEachSite)
            h += cos_mli * cos_mlj
    #
    # return h
    return round(h, 3)
#

def getForceMatrix(struct,indxNeigsFromEachSite, dRneighborsFromEachSite, radio):
    num_sites = struct.num_sites
    dim = 3 * num_sites

    # indxNeigsFromEachSite = getIndxNeigsFromEachSite(struct, radio)

    forceMatrix =  [ [ 0.5 * forceMatrix_st(s, t, num_sites, indxNeigsFromEachSite, dRneighborsFromEachSite)\
                      for s in range(dim) ]\
                      for t in range(dim) ]


    return forceMatrix # * 0.5
#
#
# s = 0
# t = 1
# indxP(s), indxI(s), indxP(t), indxI(t)
#
#
#
# radio = 1.1 # cubic
# # radio = 0.8 #fcc1.0
#
# fcc = getFccEquilibrium()
# fcc
# len(fcc.get_neighbors(fcc[0], radio, include_index=True))
# # len(fcc.get_neighbors(fcc[0], 1.1, include_index=True)) #==6 for cubic simple and alat = 1
# # len(fcc.get_neighbors(fcc[0], 0.8, include_index=True)) #==12 for fcc and alat = 1
# #


# dRneighborsFromEachSite = getDRneighborsFromEachSite(fcc, radio)
# indxNeigsFromEachSite = getIndxNeigsFromEachSite(fcc, radio)
# H = getForceMatrix(fcc,indxNeigsFromEachSite, dRneighborsFromEachSite, radio)
# print(np.matrix(H))

#
# forceMatrix_st(0, 1, fcc.num_sites, indxNeigsFromEachSite, dRneighborsFromEachSite)
#
# m = 0
# lth_neighborOfm=0
# getCos(m, lth_neighborOfm, 0, 0, indxNeigsFromEachSite, dRneighborsFromEachSite)
# getCos(m, lth_neighborOfm, 0, 1, indxNeigsFromEachSite, dRneighborsFromEachSite)
#
# import math
# from scipy import linalg as LA
#
# # For a complex Hermitian or real symmetric matrix: eigvalsh
# # getForceMatrix() returns a numpy array... OK
# # ε = LA.eigvalsh( getForceMatrix() ) # ε is a numpy array
# ε = LA.eigvalsh( H )
#
# ε
#
# H
#
#

#######
# def getEnergyOriginal(struct, radio, neighborsIndxsOfAllSites): #radio=1.5
# energy = 0
# listInterac = []
#
# # use this when atoms have departed and others have came:
# if (len(neighborsIndxsOfAllSites) == 0):
#     neighborsIndxsOfAllSites = getNeighborsIndxsOfAllSites(struct, radio)
# #
#
# for i in range(struct.num_sites):
#     site = struct[i]
#     jNeighbors = neighborsIndxsOfAllSites[i]
#
#     # sum up the interactions:
#     for j in jNeighbors:
#         # avoid repeated interactions:
#         if not ( ([i, j] in listInterac) or ([j, i] in listInterac) ):
#             energy += V( site.distance(struct[j]) )
#             listInterac.append([i, j])
# #
# return energy
# #

################################################################################
#
# def derivRderivX(n, m, p, i, struct):
#     # (X_n_i - X_m_i) / Rnm
#     if (n != m): # must be different atoms!!!
#         if (p == n):
#             return (struct[n].coords[i] - struct[m].coords[i]) / struct[n].distance[struct[m]]
#             #
#         elif (p == m):
#             return -(struct[n].coords[i] - struct[m].coords[i]) / Rnm
#     #
#     return 0
# #

# def getHnmpi(n, m, p, i, struct):
#
#
#     siteN = struct[n]
#     siteM = struct[m]
#     rN = siteN.coords
#     rM = siteM.coords
#     Xp_i = rN[p]
#
#     Xv_i = rM[n]
#     Xv_i = rM[]
#
#     Rnm = siteN.coords - siteM.coords
#     Rnm[i]
#
#     Xp_i = struct[n].p

#
#
# tupleNeigh = struct.get_neighbors(struct[i], radio, include_index=True)
#
#
# fcc[0].get_neighbors
#
# def cosine(i, j):
#     return (posNeigh[p][l+1][i] + dR[p][l][i] - x[p][i]) / r0
#
# def cosine(l, p, i, r0, posNeigh, dR, x):
#     # L = l + 1 # index neighbors begins in zero
#     # cosine = (posNeigh[p][l+1][i] + dR[p][l][i] - x[p][i]) / r0
#     # return (cosine / 6) / 2  # <<<<< 6 is the number of neighbors in 3D
#
#     # absdR = abs(dR[p][l][i])
#     # if ( ( 0.8 < absdR) and (absdR < 1.2) ):
#     #     return 0
#     # else:
#     #     return (posNeigh[p][l+1][i] + dR[p][l][i] - x[p][i]) / r0
#     #
#
#
#     return (posNeigh[p][l+1][i] + dR[p][l][i] - x[p][i]) / r0
#     ## c[l_, p_, i_, r0_] := ( (posNeigh[[p, l, i]] + dR[[p, l, i]])  -   x[[p, i]])/r0;
#
# def h1(q, j, p, i,posNeigh, dR, x, nNeighbors, rm):
#     if p == q:
#         sum = 0
#         for l in range(nNeighbors):
#             sum += cosine(l, p, j, rm, posNeigh, dR, x) * cosine(l, p, i, rm, posNeigh, dR, x)
#         #
#         return sum
#     #
#     return 0
# #
#
# def h2(q, j, p, i,posNeigh, dR, x, nNeighbors, rm, areNeb, neb):
#     if ((p != q) and areNeb[p][q]):
#         sum = 0
#         for l in range(nNeighbors):
#             L = l + 1 # index neighbors begins in zero
#             sum += -cosine(l, p, j, rm,posNeigh, dR, x) * cosine(l, p, i, rm,posNeigh, dR, x) * kroneckerDelta(neb[p][L], q)
#         #
#         return sum
#     #
#     return 0
# #
#
# def hSum(q, j, p, i,posNeigh, dR, x, nNeighbors, rm, areNeb, neb):
#     return h1(q, j, p, i,posNeigh, dR, x, nNeighbors, rm) + h2(q, j, p, i,posNeigh, dR, x, nNeighbors, rm, areNeb, neb)
# #
#
# def indxP(s):
#     return s // 3
#
# def indxI(s):
#     return s % 3
#
# def force(s, t,posNeigh, dR, x, nNeighbors, rm, areNeb, neb):
#     return hSum(indxP(s), indxI(s), indxP(t), indxI(t),posNeigh, dR, x, nNeighbors, rm, areNeb, neb)
# #
#
# def getMatrixForce():
#     fact = 72 #72 / 2
#     fact = fact * 2
#     forceMatrix = [  [fact * force(s, t, posNeigh, dR, Xeq, nNeighbors, rm, areNeb, neb) for s in range(3 * nAtoms)] for t in range(3 * nAtoms)]
#     return forceMatrix
# #
#
#
# #
