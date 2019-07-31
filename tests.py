# tests

################################################################################
# python3
import os
print(os.getcwd())
os.chdir('/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/')
print(os.getcwd())

# import graphlab
################################################################################

import mcFunctions as mc
import harmonic as ha
import copy

# def getNorm(w):
#     norm = 0
#     for i in range(len(w)):
#         norm += w[i] ** 2
#     #
#     return norm ** 0.5
# #
# def getDiff(X, Xeq):
#     dw = []
#     for j in range(len(X)):
#         # dX.append(X[i] - Xeq[i])
#         dw.append(X[j][0] - Xeq[j][0])
#         dw.append(X[j][1] - Xeq[j][1])
#         dw.append(X[j][2] - Xeq[j][2])
#     #
#     return dw
# #
# def multiply(dw, forceMatrix):
#     deltaEharmonic = 0
#     nAtoms = 27
#     dimension = 3 * nAtoms
#     for i in range(dimension):
#         for j in range(dimension):
#             deltaEharmonic += dw[i] * forceMatrix[i][j] * dw[j]
#     #
#     return deltaEharmonic


def move(i, X, Xeq, forceMatrix, L):
    X = mc.getNewConfig(i, X, L)
    rHyper, deltaEharmonic = ha.getHarmonicEnergy(X, Xeq, forceMatrix)
    return rHyper, deltaEharmonic




Xeq, a1, a2, a3, Emin, forceMatrix = ha.getConfigEquilibrium()
X = copy.deepcopy(Xeq)

Emin
E = [-162.0, -161.9, -161.8974, -161.8948]

len(forceMatrix)
id(Xeq)
id(X)

L = 0.01
random()
#######
iMax = 1000
xsol = []
ysol = []

from matplotlib.pylab import *  # for random()
import copy

w0 = 0
e0 = Emin
X0 = copy.deepcopy(X)


# n = len(forceMatrix)
# forceMatrix = [ [ int(i == j) for i in range(n)] for j in range(n)]
m=3
E[m - 1]
E[m - 2]

for i in range(iMax):
    X = copy.deepcopy(X0)

    X = mc.getNewConfig(i, X, L)
    w, e = move(i, X, Xeq, forceMatrix, L)
    e = Emin + e
    P_old2new = mc.getProbTransition(e, e0, E[m - 1], E[m - 2])
    # print(e,e0,P_old2new)
    if ( P_old2new >= random() ):
        X0 = copy.deepcopy(X)
        e0 = e
        w0 = w

    # X = mc.getNewConfig(i, X, L)
    # rHyper, deltaEharmonic = ha.getHarmonicEnergy(X, Xeq, forceMatrix)
    xsol.append( w0 )
    ysol.append( e0 - Emin )
#
from matplotlib import pyplot as plt
plt.plot(xsol, ysol, 'o')

len(ysol)







#
