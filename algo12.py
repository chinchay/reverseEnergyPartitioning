################################################################################
# python2
import os
print(os.getcwd())
os.chdir('/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/')
print(os.getcwd())
path = os.getcwd()
# import graphlab
################################################################################

import functions

import potential
from potential import VLJ as V

import numpy as np


import configuration as cf
import harmonic as ha


#########################################################333
import random
import numpy as np
from scipy.integrate import simps # Integrating using Samples
from math import log, exp
import mcFunctions5 as mc
import calc as ca

from ase.io import read # to read POSCAR for ASE

# listOfPositions = cf.getOctahedral()
# Emin = cf.getLJenergy(listOfPositions)
# forceMatrix = cf.getForceMatrix()
# use_mtp = True
# use_aseLJ = False

# options for pot: 'aseLJ' or 'mtp'
# pot = 'aseLJ'
pot = 'mlip_test.ini'
# pot = 'mlip_LJ'
################################################################################
# read from POSCAR
fileName = 'POSCAR_octahedral' #POSCAR_octahedral
cfg = read(fileName)
# struct MUST have cell=[l,l,l], and pbc=True, or it won't work with LAMMPS
cfg
################################################################################
# struct
N = cfg.get_number_of_atoms()
V = cfg.get_volume()
# listOfPositions = cf.getOctahedral()
# listOfPositions
# N = len(listOfPositions)
# V = (15)**3 # cube of 15Angstrom as in the pymatnest tutorial
# isForceMatrixRequired = False #True
# Emin, eigenvalues = cf.getLJeigenvalues(listOfPositions, epsilon, sigma, rcut, getForceMatrix)
# Emin
# eigenvalues
# epsilon = 0.1
# sigma   = 2.5
# rcut    = 7.50
# aCell   = 15.0

model = cf.Model(path, pot, cfg)

# model.getEnergy(struct)


Emin, eigen = model.getEnergyAndEigen(cfg)
Emin
eigen


Emin
eigen
######

# X = struct.get_positions()
# Emin, eigen = cf.getEnergyAndEigen(X, getForceMatrix)

# if use_mtp:
#     Emin, eigenvalues = cf.getLJeigenvalues2(listOfPositions, epsilon, sigma, rcut, getForceMatrix, aCell)
# elif use_aseLJ:
#     Emin, eigenvalues = cf.getLJeigenvalues(listOfPositions, epsilon, sigma, rcut, getForceMatrix)
#
#
# Emin
# eigenvalues

# Emin = cf.getLJenergy(listOfPositions, sigma, epsilon)
# Emin

maxDeltaE = 0.5 * sum(eigenvalues) * ( 0.01 ** 2)
maxDeltaE

deltaE = maxDeltaE #/ 100
E1 = Emin + deltaE
E2, I1, I2, errorPercent = ca.getEnext(eigenvalues, N, V, Emin, E1)
E2, I1, I2, errorPercent


# Plot
from matplotlib import pyplot as plt
intervalos = 50
x  = np.array( [ (i * deltaE / intervalos) for i in range(intervalos + 8) ] )
y  = cf.harmonicDOS(x, eigenvalues, N, V) # x departs from Emin -> interpreted as dE
y1 = cf.harmonicDOS(deltaE, eigenvalues, N, V)
y2 = cf.harmonicDOS(E2-Emin, eigenvalues, N, V)


x = x + Emin
plt.plot(x, y, '-',   [Emin, E1, E2], [0,y1,y2], 'o')


x2 = np.array( [ x[i] for i in range( 20, len(x) )] )
y2 = np.array( [ log(y[i]) for i in range( 20, len(y) )] )
len(x2)
len(y2)
plt.plot(x2, y2, '-',)


################################################################################
log2 = log(2)  #math.log(2.0)




################################################################################

#because we are dealing with intervals:
# E_m values:
E            = [Emin, E1, E2] # = [E1, E2]
iDOS         = [0,     1,  1] # IntegralDosFromEminToE1 = I1, IntegralDosFromE1ToE2 = I2, and I1=I2
# log_idos     = [0, 0] # just the Ln of idos
# log_sum_idos = log(1 + 1) # = log( sum(exp(log_idos))) = log(1 + 1)
lAlpha       = [0,     0,  1]     # alpha1 = I2 / I1
# [α1, α2, ..., αm, ...] #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
m = 2
################################################################################





# Because python lists begin with index in zero, then `m` will begin
# in zero, instead of 1.

# log_idos = [0, 0] # = [log_idos_m=0, log_idos_m=1] #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# log_idos = [0, 0] # = [log_idos_m=0, log_idos_m=1]
# log_sum_idos = log(1 + 1) # = log( sum(exp(log_idos))) = log(1 + 1)

# α_m=0 = 1: impossed
# α_m=1 = 1: impossed
# α_m=2 = 1: because log_idos[0] = log_idos[1]
# lAlpha = [1, 1, 1] # [α0, α1, α2, ..., αm, ...]
# lAlpha = [1, 1, 1] # [α0, α1, α2, ..., αm, ...] #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# m = 2  # python begins in 0. << the number of the next subdivision (instead of m=3 as in the paper)
#m = 3  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


################################################################################
# E_m values:
# E = [Emin, E1, E2] # = [E0, E1, E2]
E

#test
E3temp, I1, I2, errorPercent = ca.getEnext(eigenvalues, N, V, Emin, E2)
E3temp, I1, I2, errorPercent
[Emin - Emin, E1 - Emin, E2 - Emin, mc.getEm(lAlpha[m], E[m], E[m - 1]) - Emin]
[Emin - Emin, E1 - Emin, E2 - Emin, E3temp - Emin]
################################################################################


################################################################################
################################################################################
continuar = True
iMax   = 200 #30 #10 #100 #15 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
nSteps =  40000 #20000 #5000 #20000 #1000 #5000 # 500 #50#  = nBarridos * nAtomsInConfig
alphaMin = 0.01#0.5 # the same as the paper. "At this point, the newly guessed Em is
               # deleted and then a set of random (all moves are accepted) MC
               # moves are performed... more info about how to do the last part.


# Main Reverse Energy Partitioning loop:
# L = 0.0001
L = 0.06 #0.01 #~amplitud of random walk # usado para el LJ de ASE<<<<<<<<<<<<<<<<<<<<<
# L = 0.14
# L = 0.12 #0.01 #~amplitud of random walk
# L = 0.28
# L=0.06
# L =0.008
# L=0.004
# L = 0.05
# L = 0.12
# L0 = L

i = 0
# nSteps0 = nSteps


import copy
e = Emin
Xeq = listOfPositions
X = copy.deepcopy(Xeq)
L0 = L

################################################################################



# X = get
# e = mc.getEnergyConfig(X)
from matplotlib.pylab import *  # for random()
def getRandomWalk(L): # L: amplitud
    return L * (random() - 0.5) # `0.5`:for centering around 0.
#
import copy
import os
import vectormath as vmath
from vectormath import Vector3 as Vector3
from functions import belongs




# MC sampling. Collect quantities for the next subdivision m+1.
# Xeq = copy.deepcopy(X)

# def unitVect(i):
#     if i == 0:
#         return  Vector3(1,0,0)
#     elif i == 1:
#         return  Vector3(0,1,0)
#     elif i == 2:
#         return  Vector3(0,0,1)
# #
def unitVect2(i):
    if i == 0:
        return  Vector3(1,0,0)
    elif i == 1:
        return  Vector3(0,1,0)
    elif i == 2:
        return  Vector3(0,0,1)
    elif i == 3:
        return  Vector3(-1,0,0)
    elif i == 4:
        return  Vector3(0,-1,0)
    elif i == 5:
        return  Vector3(0,0,-1)
#
def s(n,L):
    DimTotal = n * 3
    r = random() * DimTotal
    iDim = int(r )
    iAtom = iDim // 3
    dim =  iDim % 3
    l = ( r - ( (iAtom * 3) + dim) ) - 0.5
    return iAtom, dim, l*L
#

# Static variable, for counting:
def counter():
    if 'cnt' not in counter.__dict__:
        counter.cnt = 0
    counter.cnt += 1
    return counter.cnt
#

def draw(X, e):
    c = counter()
    # if (c % 5 == 0):
    # if (c % 100 == 0):
    if (c % 1000 == 0):
        nAtoms = len(X)
        out = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/moves." + str(c) + ".xyz"
        pos = str(nAtoms) + "\n\n"
        for i in range(len(X)):
            pos += "H   " + str(X[i][0]) + "   " + str(X[i][1]) + "   " + str(X[i][2]) + "\n"
        #
        pos += "\n"

        f = open(out, "w")
        f.write(pos)
        f.close()

        out = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/configs.txt"
        string = "moves." + str(c) + ".xyz" + " , " + str(e) + "\n"
        f = open(out, "a")
        f.write(string)
        f.close()
#

def maintainCenterOfMass(X):
    r0  = Vector3(0,0,0)
    rCM = Vector3(0,0,0)
    n = len(X)
    for iAtom in range(n):
        rCM += X[iAtom]
    #
    rCM = rCM / n # <<<<<<<<<<<<<, include mass of atoms later <<<<<<<<<<<<<<<<<<<<<<<
    for iAtom in range(n):
        X[iAtom] -= rCM
    #
    return X
#
# def maintainCenterOfMass3(cfg):
#     return cfg.translate(-cfg.get_center_of_mass())
    #
#


def translate(Xold, L):
    Xnew = copy.deepcopy(Xold)
    iAtom, idim, l = s(6, L)
    d = unitVect2(idim) * L * (random() - 0.5)
    Xnew[iAtom] += d
    return Xnew
#

# def getInitialStruct():
#     import vectormath as vmath
#     from vectormath import Vector3 as Vector3
#     l = 15.0
#     a1 = Vector3(5.35391903, 7.34897441, 9.20345327)
#     a2 = Vector3(4.83487693, 9.69665956, 7.78048212)
#     a3 = Vector3(3.48505708, 8.10763956, 11.13648625)
#     a4 = Vector3(2.73813043, 7.97569725, 8.44800613)
#     a5 = Vector3(5.58070460, 9.82836905, 10.46934040)
#     a6 = Vector3(2.96585081, 10.45497144, 9.71365159)
#     positions = [a1, a2, a3, a4, a5, a6]
#     #
#     n1 = len(listOfPositions)
#     from ase import Atom, Atoms
#     chem = 'H' + str(n1)
#     cfg = Atoms( chem, positions=positions , cell=(l,l,l), pbc=True )
#     #
#     return positions, cfg
# #
#
# def updateCfgPositions(cfg, X):
#     from ase import Atom, Atoms
#     positions = [ (X[i][0], X[i][1], X[i][2]) for i in range(X)]
#     cfg.set_positions( positions )
#     return cfg
# #


def swap(X, neighbors, S): #S:chemical symbols list
    # cfg = updateCfgPositions(cfg, X)
    Xnew = copy.deepcopy(X)
    Snew = copy.deepcopy(S)
    #
    i  = random.randint(0, nAtoms)
    j  = random.randint(0, len(neighbors))
    #
    # assert( cfg.get_chemical_symbols()[i] != cfg.get_chemical_symbols()[j] )
    assert( S[i] != S[j] )

    xi = X[i].copy()
    xj = X[j].copy()
    si = S[i].copy()
    sj = S[j].copy()
    Xnew[i] = xj
    Xnew[j] = xi
    Snew[i] = sj
    Snew[j] = si
    # Cfgnew = updateCfgPositions(cfg, X)
    return Xnew, Snew
#

def getMove(nAtoms, X, canSwap, neighbors, S): # positionsVec, S:chemSymbols
    import random
    # i = random.randint(0, nAtoms)
    # we can translational move Xi, or swap with a neighbor
    if canSwap:
        if random() < 0.5: # random translational move:
            Xnew = translate(X, L)
            Snew = S # no problem not doing deepcopy(), since it is the same, no swap.
        else: # random swapping move:
            Xnew, Snew = swap(X, neighbors, S)
    else:
        Xnew = translate(X, L)
    #
    return Xnew, Snew
#

canSwap = True
N3 = 3 * N
import matscipy.neighbours
r_cut = 2.8

def getRandomNeighbor(cfg, i):
    Z  = cfg.get_atomic_numbers()
    (i_list, j_list) = matscipy.neighbours.neighbour_list('ij', cfg, r_cut) # r_cut is glabal variable. Should it be local?
    indx_list = np.where( ( i_list == i ) & ( Z[ i_list ] != Z[ j_list ] ) )[0]
    n = len(indx_list)
    if n >= 1:
        return j_list[ indx_list[ random.randint(0, n - 1) ]  ]
    #
    return -1
#
def translate3(cfg, L):
    i = random.randint(0, N3)
    cfgTemp = copy.deepcopy(cfg)
    cfgTemp.positions[ i // 3 ][ i % 3 ] += L * (random() - 0.5)
    return cfgTemp
#
def swap3(cfg):
    i = random.randint(0, N)
    j = getRandomNeighbor(cfg, i)
    cfgTemp = cfg
    if j != -1:
        cfgTemp = copy.deepcopy(cfg)
        cfgTemp.positions[i] = copy.deepcopy( cfg.positions[j] )
        cfgTemp.positions[j] = copy.deepcopy( cfg.positions[i] )
        return cfgTemp
    #
    return cfgTemp
#


def getMove3(cfg): # positionsVec, S:chemSymbols
    # i = random.randint(0, nAtoms)
    # we can translational move Xi, or swap with a neighbor
    if canSwap:
        if random() < 0.5: # random translational move:
            cfgTemp = translate3(cfg, L)
        else: # random swapping move:
            cfgTemp = swap3(cfg)
    else:
        cfgTemp = translate3(cfg, L)
    #
    return cfgTemp
#

# nSteps
# N


def move3(cfg, nSteps, Emin, Em, Emplus1, L):
    hasMoved = False
    acceptances = 0
    [I_EminToEm, I_EmToEmplus1, I_Emplus1ToInfty] = [0, 0, 0]
    Eold = model.getEnergy(cfg)
    for i in range(nSteps):
        # iAtom, idim, l = s(6, L)
        if hasMoved:
            hasMoved = False
            Eold = e
        #
        cfgTemp = getMove3(cfg)
        e = model.getEnergy(cfgTemp)
        if (e <= Emplus1):
            hasMoved = True
            cfg = copy.deepcopy(cfgTemp)
            acceptances += 1 # accepted move.
            draw(cfg.positions, e)
            cfg.translate( -cfg.get_center_of_mass() )

        else:
            e = Eold
        #
        if belongs(e, Emin, Em): # (Emin <= e <= E[m - 1]): #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            I_EminToEm += 1  # un gol más para ∫_Emin^Em-1 Ω(E)dE
        elif belongs(e, Em, Emplus1): # (E[m - 1] <= e <= Em):
            I_EmToEmplus1 += 1 # un gol más para ∫_Em-1^Em Ω(E)dE
        elif Emplus1 < e:
            I_Emplus1ToInfty += 1
        #
    #
    ratioAcceptances = float(acceptances) / nSteps #`float()` to ensure the division will be a float (python2 resembles Fortran in divisions)

    if I_EminToEm == 0:
        I_EminToEm = 1
    alpha = I_EmToEmplus1 / float(I_EminToEm) # `float()` to ensure the division will be a float (python2 resembles Fortran in divisions)
    print("alpha = ", float("{0:.3f}".format(alpha)))
    print(float("{0:.3f}".format(ratioAcceptances)), I_EminToEm, I_EmToEmplus1, I_Emplus1ToInfty, "...", L)
    return alpha, ratioAcceptances
#



def move2(X, sigma, epsilon, rcut, nSteps, Emin, Em, Emplus1, L):
    hasMoved = False
    acceptances = 0
    [I_EminToEm, I_EmToEmplus1, I_Emplus1ToInfty] = [0, 0, 0]
    Eold = cf.getLJenergy(X, sigma, epsilon)
    for i in range(nSteps):
        # iAtom, idim, l = s(6, L)
        if hasMoved:
            hasMoved = False
            Eold = e
        #
        # d = unitVect2(idim) * L * (random() - 0.5)
        ########################################################################
        Xt, St = getMove(nAtoms, X, canSwap, neighbors, S) #move includes translation or swap
        # Xt[iAtom] += d
        ########################################################################

        ########################################################################
        #
        getForceMatrix = False
        #
        if use_mtp:
            e, eigenvalues = cf.getLJeigenvalues2B(Xt, St, epsilon, sigma, rcut, getForceMatrix, aCell)
        elif use_aseLJ:
            e, eigenvalues = cf.getLJeigenvaluesB(Xt, St, epsilon, sigma, rcut, getForceMatrix)
        #
        if (e <= Emplus1):
            hasMoved = True
            # X = copy.deepcopy(Xt)

            # X[iAtom] += d
            # X = copy.deepcopy(Xt)
            # S = copy.deepcopy(St)
            X = Xt # ???????????????????????????????????????????????????????????
            S = St

            acceptances += 1 # accepted move.
            draw(X, e)
            X = maintainCenterOfMass(X)
        else:
            e = Eold
        #
        #
        if belongs(e, Emin, Em): # (Emin <= e <= E[m - 1]): #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            I_EminToEm += 1  # un gol más para ∫_Emin^Em-1 Ω(E)dE
        elif belongs(e, Em, Emplus1): # (E[m - 1] <= e <= Em):
            I_EmToEmplus1 += 1 # un gol más para ∫_Em-1^Em Ω(E)dE
        elif Emplus1 < e:
            I_Emplus1ToInfty += 1
        #
    #
    ratioAcceptances = float(acceptances) / nSteps #`float()` to ensure the division will be a float (python2 resembles Fortran in divisions)

    if I_EminToEm == 0:
        I_EminToEm = 1
    alpha = I_EmToEmplus1 / float(I_EminToEm) # `float()` to ensure the division will be a float (python2 resembles Fortran in divisions)
    print("alpha = ", float("{0:.3f}".format(alpha)))
    print(float("{0:.3f}".format(ratioAcceptances)), I_EminToEm, I_EmToEmplus1, I_Emplus1ToInfty, "...", L)
    return alpha, ratioAcceptances
#

def move(X, sigma, epsilon, rcut, nSteps, Emin, Em, Emplus1, L):
    hasMoved = False
    acceptances = 0
    [I_EminToEm, I_EmToEmplus1, I_Emplus1ToInfty] = [0, 0, 0]
    Eold = cf.getLJenergy(X, sigma, epsilon)
    # acceptances
    # I_EminToEm
    # I_EmToEmplus1
    for i in range(nSteps):
        # iAtom = i % len(X)
        # iAtom = int(random() * len(X))
        # idim  = int(random() * 6)

        # get a configuration e,X after randomly move X:
        # Eold = getEnergyConfig(X)
        # assert( Eo == getEnergyConfig(X) )

        iAtom, idim, l = s(6, L)

        if hasMoved:
            hasMoved = False
            Eold = e
        #
        # Xt = copy.deepcopy(X)

        # d = Vector3(getRandomWalk(L), getRandomWalk(L), getRandomWalk(L))
        # d = unitVect(idim) * getRandomWalk(L)
        d = unitVect2(idim) * L * (random() - 0.5)
        # d = unitVect2(idim) * l


        ########################################################################
        Xt[iAtom] += d
        ########################################################################

        ########################################################################
        #
        getForceMatrix = False
        #
        if use_mtp:
            e, eigenvalues = cf.getLJeigenvalues2(Xt, epsilon, sigma, rcut, getForceMatrix, aCell)
        elif use_aseLJ:
            e, eigenvalues = cf.getLJeigenvalues(Xt, epsilon, sigma, rcut, getForceMatrix)
        #
        # e  = cf.getLJenergy(Xt, sigma, epsilon)


        # P_old2new = mc.getProbTransition(e, Eold, Em, EmMinus1)
        # P_old2new = mc.getProbTransition(e, Eold, Emplus1 + 0.0004, Em)
        # P_old2new = mc.getProbTransition(e, Eold, Emplus1, Em)

        # os.system("echo "+ str(e) + " >> out2.txt")

        # if (e <= Emplus1 + 0.001):
        # if (e <= Emplus1 + 0.0006):
        # if (e <= Emplus1 + 0.0002):
        if (e <= Emplus1):
        # if ( P_old2new >= random() ):
            # os.system("echo "+ str(e) + " >> out.txt")
            hasMoved = True
            # X = copy.deepcopy(Xt)
            X[iAtom] += d
            acceptances += 1 # accepted move.
            draw(X, e)
            X = maintainCenterOfMass(X)
        else:
            e = Eold
        #

        # # it has no sense to add the same configuration (not moved) to Ω:
        # if hasMoved:
        #     acceptances += 1 # accepted move.
        # #
        if belongs(e, Emin, Em): # (Emin <= e <= E[m - 1]): #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            I_EminToEm += 1  # un gol más para ∫_Emin^Em-1 Ω(E)dE
        elif belongs(e, Em, Emplus1): # (E[m - 1] <= e <= Em):
            I_EmToEmplus1 += 1 # un gol más para ∫_Em-1^Em Ω(E)dE
        elif Emplus1 < e:
            I_Emplus1ToInfty += 1
        #
    #
    ratioAcceptances = float(acceptances) / nSteps #`float()` to ensure the division will be a float (python2 resembles Fortran in divisions)

    if I_EminToEm == 0:
        I_EminToEm = 1
    alpha = I_EmToEmplus1 / float(I_EminToEm) # `float()` to ensure the division will be a float (python2 resembles Fortran in divisions)
    print("alpha = ", float("{0:.3f}".format(alpha)))
    print(float("{0:.3f}".format(ratioAcceptances)), I_EminToEm, I_EmToEmplus1, I_Emplus1ToInfty, "...", L)
    # acceptances
    # I_EminToEm
    # I_EmToEmplus1
    ################################################################################

    # log_idos.append( log_sum_idos + log(alpha) ) # =log_idos[m]
    #
    # # calculate log_sum_idos:
    # # Log of the sum of integrals for each interval, m=1, m=2, m=3,...
    # DOSintegrals = [ exp(log_int) for log_int in log_idos ]
    # log_sum_idos = log(sum(DOSintegrals))
    #
    # return alpha, log_idos, log_sum_idos, L
    # return alpha
    return alpha, ratioAcceptances
#

def goodMove(Xo, sigma, epsilon, rcut, nSteps, Emin, Em, Emplus1, L):
    ratioMax = 0.7 #0.35
    ratioMin = 0.1 #0.25
    good = False
    count = 0
    while ((not good) and count <= 10):
        count += 1
        alpha, ratioAcceptances = move(Xo, sigma, epsilon, rcut, nSteps, Emin, Em, Emplus1, L)
        if ratioAcceptances > ratioMax:
            L = 1.1 * L
            print("...... increasing L = ", L, ratioAcceptances)
            X = copy.deepcopy(Xo)
            # return 0
        elif ratioAcceptances < ratioMin:
            L = 0.9 * L
            print("...... decreasing L = ", L, ratioAcceptances)
            X = copy.deepcopy(Xo)
            # return 0
        else:
            print("good.")
            good = True
    #
    # log_idos.append( log_sum_idos + log(alpha) ) # =log_idos[m]
    # DOSintegrals = [ exp(log_int) for log_int in log_idos ]
    # log_sum_idos = log(sum(DOSintegrals))

    # if good:
        # log_idos.append( log_sum_idos + log(alpha) ) # =log_idos[m]

        # calculate log_sum_idos:
        # Log of the sum of integrals for each interval, m=1, m=2, m=3,...
        # DOSintegrals = [ exp(log_int) for log_int in log_idos ]
        # log_sum_idos = log(sum(DOSintegrals))

    return alpha, ratioAcceptances, L #, log_idos, log_sum_idos, L
#

#=======================================================

# EmMinus1 = E[m - 1]
# Emplus1 = mc.getEm(lAlpha[m], E[m], E[m - 1])
# Emplus1
# E.append(Emplus1)
# E
# [E[i]-Emin for i in range(len(E))]
#
# X = copy.deepcopy(Xeq)
# Em = E[m]
# Em
# m
# # os.system("rm out.txt")
# # os.system("rm out2.txt")
# # os.system("touch out.txt")
# # os.system("touch out2.txt")
#
# [E[i]-Emin for i in range(len(E))]
#

################################################################################
# Eo  = getEnergyConfig(X)
# Xo = copy.deepcopy(X)





#
#
# lCfgs, alpha, log_idos, log_sum_idos, L =\
#                     mc.randomMCmoves(Emin, E[m - 1], E[m], E[m + 1]+0.00,\
#                                   nSteps, e, Xcopy, log_idos, log_sum_idos, L)





# import os
# print(os.getcwd())
# os.chdir('/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/')
# print(os.getcwd())
# import numpy as np
# from matplotlib import pyplot
# columns = np.loadtxt("out.txt", unpack=True)
# columns = columns - Emin
# bin_contents, bin_ends, patches = pyplot.hist( columns, bins = 100  , range     = (min(columns), max(columns))       )
#
#
# columns = np.loadtxt("out2.txt", unpack=True)
# columns = columns - Emin
# bin_contents, bin_ends, patches = pyplot.hist( columns, bins = 400  , range     = (min(columns), max(columns))       )
# bin_width     = (bin_ends[1]-bin_ends[0])
# bin_widths    = len(bin_contents)*[bin_width]
# bin_positions = bin_ends[0:-1]+bin_width/2.0

################################################################################
# iii
Xo = copy.deepcopy(Xeq)
while ( continuar and (i < iMax) ): # iMax allows to force stopping.
    i += 1

    EmMinus1 = E[m - 1]
    Em       = E[m]
    Emplus1  = mc.getEm(lAlpha[m], Em, EmMinus1)
    E.append(Emplus1)


    # [ Emin-Emin, E[m - 1]-Emin, E[m]-Emin, E[m + 1] ]

    # X = get
    # e = mc.getEnergyConfig(X)

    # MC sampling. Collect quantities for the next subdivision m+1.
    alpha, ratioAcceptances, L = goodMove(Xo, sigma, epsilon, rcut, nSteps, Emin, Em, Emplus1, L)

    # Xcopy = copy.deepcopy(X)
    # lCfgs, alpha, log_idos, log_sum_idos, L =\
    #                     mc.randomMCmoves(Emin, E[m - 1], E[m], E[m + 1],\
    #                                   nSteps, e, Xcopy, log_idos, log_sum_idos, L)

    #

    # if ((alpha < 0.5) and (i < 30)):
    #     print("because alpha is < 0.5, we will repeat the step.")
    #     print(E[m])
    #     E[m] = E[m] + 0.1
    #     # nn = len(E) - 1
    #     # Em = E[nn]
    #     # print(Em)
    #     # # Em += 0.1 #((E[m - 1] - E[m - 2]) * 10)
    #     # del E[nn]
    #     # E.append(Em)
    #     m -= 1
    # #
    #



    # L = L0
    # nSteps = nSteps0
    lAlpha.append(alpha) # alpha = α_m
    print(".............ratioAcceptances, alpha, Energy:", m, ratioAcceptances, Emplus1, "alpha = ", alpha )

    # select a starting configuration for the next subdivision:
    # assert( len(lCfgs) != 0 )
    # nthMin = 1
    # e, X = mc.select_start_config(lCfgs, E[m - 1], E[m], nthMin)

    #---------------
    m += 1 # update the number of energy subdivisions m:
    #---------------

    # adjust MC step size for the nex subdivision:
    #nSteps = mc.adjust_step_size() # step_size (=L) or nSteps ???

    # finish as α < α_min
    if (alpha < alphaMin):
        continuar = False # it will stop the main loop.
        m -= 1 # to redo the m-th partition.

        Einfty = E[m - 1] + 10000 # Einfty -> \infty # "all moves are accepted"
        # I think E[m] is infty now... should I overwrite E[m]????? <<<<<<<<<<<<

        # do a set of random MC moves and collect ehist above and below E[m-1]:
        alpha, ratioAcceptances, L = goodMove(Xo, sigma, epsilon, nSteps, Emin, Em, Einfty, L)

        # lCfgs, alpha, log_idos, log_sum_idos, L =\
        #             mc.randomMCmoves(Emin, E[m - 2], E[m - 1], Einfty,\
        #                           nSteps, e, X, log_idos, log_sum_idos, L)

        lAlpha[m] = alpha # <- overwrite since α_m was recalculated.
    #
#
L
Esalvado = E
lAlphaSalvado = lAlpha
len(E)
len(lAlpha)
# m #97
# m=len(E)-1

iDOS = [1, 1] # IntegralDosFromEminToE1 = I1, IntegralDosFromE1ToE2 = I2, and I1=I2
iDOS_ideal = [1, 1]
ints = [ i for i in range(len(lAlpha) - 3) ]
for i in range(len(ints)):
    iDOS.append( lAlpha[i + 3] * sum(iDOS) )
    iDOS_ideal.append( 1 * sum(iDOS) )
#
ints = [ i for i in range(len(lAlpha) - 1) ]
len(ints)
len(iDOS)
plt.plot( ints, [log(iDOS[i]) for i in range(len(ints))], 'o', ints, [log(iDOS_ideal[i]) for i in range(len(ints))], 'o' )

lAlpha[0 + 3]
lAlpha


################################################################################
# normalize Ω (taking advantage of logarith property when substracting):
# log_idos = [log_idos[i] - log_sum_idos for i in range(len(log_idos)) ]

# len(E)
# len(iDOS)
# len(lAlpha)

# le = min(len(E), len(lAlpha))

for m in range(len(lAlpha)):
    print( E[m], lAlpha[m] )
#

# l = len(E)
for i in range(len(E)):
    if (i + 1 < len(E)): print(E[i + 1] - E[i])

# ll = le
plt.plot([lAlpha[i + 3] for i in range(len(lAlpha) - 3)], '+')

len(E)
len(Eshifted)
len(lAlpha)
Eshifted = [ E[i] - Emin for i in range(len(E)) ]
plt.plot([Eshifted[i+3] for i in range(len(Eshifted)-3)], [lAlpha[i+2] for i in range(len(lAlpha)-2)], '+')

log_Energies     = [ log( E[i+1] - Emin )    for i in range(len(iDOS)) ]
plt.plot([log_Energies[i+1] for i in range(len(log_Energies)-1)], [lAlpha[i+2] for i in range(len(lAlpha)-2)], '+')

log_Integral_DOS_EminToEi = [ log( sum(iDOS[0:i+1]) ) for i in range(len(iDOS)) ]

# iDOS_harmonic = [ cf.harmonicDOS(Eshifted[i+1], eigenvalues, N, V) ]
const = cf.Integral_harmonicDOS(1.0, eigenvalues, N, V)
const
# import imp
# imp.reload(cf)

log_Integral_harmonicDos_EminToEi = [ log( cf.Integral_harmonicDOS(Eshifted[i+1], eigenvalues, N, V) / (const**3)) for i in range(len(iDOS)) ]

plt.plot( log_Energies, log_Integral_DOS_EminToEi, 'o', log_Energies, log_Integral_harmonicDos_EminToEi , 'o')


np.polyfit(log_Energies[0:20], log_Integral_DOS_EminToEi[0:20], 1)
a, b = np.polyfit(log_Energies[0:20], log_Integral_harmonicDos_EminToEi[0:20], 1)
D = a * 2
D
N = (D + 3) / 3
N
N = (D + 6) / 3
N

np.polyfit(log_Energies[0:20], log_Integral_DOS_EminToEi[0:20], 1)
np.polyfit(log_Energies[0:20], log_Integral_harmonicDos_EminToEi[0:20], 1)
eigenvalues



# len(iDOS)
# len(E)
# ll
# E_b    = E[0:ll]
# iDOS_b = iDOS[0:ll]
# from scipy import interpolate
# f = interpolate.interp1d(E_b, iDOS_b)
# xnew = np.arange(min(E_b), max(E_b), 0.001)
# ynew = f(xnew)
# plt.plot(E_b, iDOS_b, 'o', xnew, ynew, '-')
#
# from scipy.misc import derivative
# xnew2 = xnew[1:len(xnew)]
# df = [derivative(f, x, dx=1e-6, n=1) for x in xnew2 ]
# plt.plot(xnew2, df, '-')
#
# import analyze as an
# # le = len(df)
# # le
# partitions = 1000
# Tmin = 10
# Tmax = 1000
# k_Boltzmann   = 8.617
#
# listT   = [ Tmin + (i * ( (Tmax - Tmin) / partitions )) for i in range(partitions) ]
# listQex = [ an.getQex( 1/(listT[i] * k_Boltzmann), xnew2, df) for i in range(len(listT)) ]
# plt.plot(listT, listQex, '-')
#
#
# len(E)
# len(log_idos)
# E[0]
# log_idos[0]
#
#
# dos_meanValues = [ iDOS[i] / (E[i] - E[i-1]) for i in range(ll) ]
# E_medios       = [ (E[i] + E[i-1]) / 2 for i in range(ll) ]
# beta = k_Boltzmann * T
# y1 = [ dos_meanValues[i] * exp(-E_medios[i] * beta) for i in range(ll) ]
# y2 = [ E_medios[i] * y1[i] for i in range(ll) ]
# y3 = [ (E_medios[i] ** 2) * y1[i] for i in range(ll) ]
#
#
def getDOSmeanValues(E, iDOS):
    assert( len(E) == len(iDOS) + 2 )
    return [ iDOS[i] / (E[i+1] - E[i]) for i in range(len(iDOS)) ]
#

def getEmedios(E):
    return [ (E[i+1] + E[i]) / 2 for i in range(len(E) - 1) ]

# def getPartitionFunctionZ(k_Boltzmann, T, E, iDOS): # T : temperature
#     assert( len(E) == len(iDOS) + 1 )
#     dos_meanValues = getDOSmeanValues(E, iDOS)
#     E_medios       = getEmedios(E)
#     beta = k_Boltzmann * T
#     Z = 0
#     for i in range(len(iDOS)):
#         Z += dos_meanValues[i] * exp(-E_medios[i] * beta)
#     #
#     return Z
# #
#
# def getExpectedEnergy(k_Boltzmann, T, E, iDOS): # T : temperature
#     assert( len(E) == len(iDOS) + 1 )
#     dos_meanValues = getDOSmeanValues(E, iDOS)
#     E_medios       = getEmedios(E)
#     beta = k_Boltzmann * T
#     expectedEnergy = 0
#     for i in range(len(iDOS)):
#         expectedEnergy += E_medios[i] * dos_meanValues[i] * exp(-E_medios[i] * beta)
#     #
#     return expectedEnergy / getPartitionFunctionZ(k_Boltzmann, T, E, iDOS)
# #
#
# def getExpectedEnergySquared(k_Boltzmann, T, E, iDOS): # T : temperature
#     assert( len(E) == len(iDOS) + 1 )
#     dos_meanValues = getDOSmeanValues(E, iDOS)
#     E_medios       = getEmedios(E)
#     beta = k_Boltzmann * T
#     expectE2 = 0
#     for i in range(len(iDOS)):
#         expectE2 += (E_medios[i] ** 2) * dos_meanValues[i] * exp(-E_medios[i] * beta)
#     #
#     return expectE2 / getPartitionFunctionZ(k_Boltzmann, T, E, iDOS)
# #

def getExpectedValues(k_Boltzmann, T, E, iDOS): # T : temperature
    assert( len(E) == len(iDOS) + 2 )
    dos_meanValues = getDOSmeanValues(E, iDOS)
    E_medios       = getEmedios(E)
    beta = 1 / (k_Boltzmann * T)
    expectedEnergy = 0
    expectedEnergySq = 0
    partitionFunctionZ = 0
    for i in range(len(iDOS)):
        prob = exp(-E_medios[i] * beta) * dos_meanValues[i]
        partitionFunctionZ += prob
        expectedEnergy     += prob * E_medios[i]
        expectedEnergySq   += prob * (E_medios[i] ** 2)

    #
    expectedEnergy   = expectedEnergy / partitionFunctionZ
    expectedEnergySq = expectedEnergySq / partitionFunctionZ
    return expectedEnergy, expectedEnergySq

def getSpecificHeat(k_Boltzmann, T, E, iDOS): # T : temperature
    expectedEnergy, expectedEnergySq = getExpectedValues(k_Boltzmann, T, E, iDOS)
    return (1 / (k_Boltzmann * (T ** 2))) * ( expectedEnergySq - (expectedEnergy ** 2) )
#
N
# E = E[0:98]
# E_new = [E[i] - Emin for i in range(len(E))]
# E_new
len(Eshifted)
len(iDOS)
# iDOS = iDOS[0:97]
# iDOS_harmonic = [ cf.harmonicDOS(E[i+1], eigenvalues, N, V) - cf.harmonicDOS(E[i], eigenvalues, N, V) for i in range(len(E) - 1) ]
# k_Boltzmann = 8.617/100000 # units eV/Kelvin
k_Boltzmann = 8.617e-5 # units eV/Kelvin

iDOS_harmonic = [ cf.Integral_harmonicDOS(Eshifted[i+1], eigenvalues, N, V) - cf.Integral_harmonicDOS(Eshifted[i], eigenvalues, N, V) for i in range(len(Eshifted)-2) ]
len(iDOS_harmonic)

temperatures = [i for i in range(1,1000)]
listExpectedEnergies = [getExpectedValues(k_Boltzmann, temperatures[i], Eshifted, iDOS)[0] for i in range(len(temperatures))]
listExpectedEnergiesHarmonic = [getExpectedValues(k_Boltzmann, temperatures[i], Eshifted, iDOS_harmonic)[0] for i in range(len(temperatures))]
plt.plot(temperatures, listExpectedEnergies, temperatures, listExpectedEnergiesHarmonic)

Emin
listExpectedEnergies

temperatures = [i*10 for i in range(1,1000)]
listExpectedEnergies = [getExpectedValues(k_Boltzmann, temperatures[i], Eshifted, iDOS)[0] for i in range(len(temperatures))]
listExpectedEnergiesHarmonic = [getExpectedValues(k_Boltzmann, temperatures[i], Eshifted, iDOS_harmonic)[0] for i in range(len(temperatures))]
plt.plot(temperatures, listExpectedEnergies, temperatures, listExpectedEnergiesHarmonic)

temperatures = [i/float(20) for i in range(1,1000)] #temperatures = [i/float(50) for i in range(1,1000)]
listExpectedEnergies = [getExpectedValues(k_Boltzmann, temperatures[i], Eshifted, iDOS)[0] for i in range(len(temperatures))]
listExpectedEnergiesHarmonic = [getExpectedValues(k_Boltzmann, temperatures[i], Eshifted, iDOS_harmonic)[0] for i in range(len(temperatures))]
plt.plot(temperatures, listExpectedEnergies, temperatures, listExpectedEnergiesHarmonic, [0],[0])

temperatures = [i for i in range(10,1000)] #temperatures = [i for i in range(3,1000)]
listCv = [getSpecificHeat(k_Boltzmann, temperatures[i], Eshifted, iDOS) for i in range(len(temperatures))]
listCvHarmonic = [getSpecificHeat(k_Boltzmann, temperatures[i], Eshifted, iDOS_harmonic) for i in range(len(temperatures))]
plt.plot(temperatures, listCv, temperatures, listCvHarmonic)

################################################################################








################################################################################
# you have finished. Now you have `E` (list of Em=i) and log_idos (log(DOS_i))
################################################################################
path = os.getcwd()
fileName = "out.txt"
f = open(path + "/" + fileName, "w")
for m in range(len(E)):
    f.write(E[m], log_idos[m], lAlpha[m])
#
f.close()
################################################################################
#
#
path = os.getcwd()
file2 = "file.energies"
f = open(path + "/" + file2, "w")
for m in range(len(E)):
    f.write(E[m], log_idos[m], lAlpha[m])
#
f.close()


################################################################################
import imp
imp.reload(an)
import analyze as an
le = min(len(E), len(log_idos))

partitions = 1000
Tmin = 10
Tmax = 300
k_Boltzmann   = 1

listEnergy    = [ E[i] for i in range(le)]
listDOS       = [ exp( log_idos[i] ) for i in range(le) ]
listT   = [ Tmin + (i * ( (Tmax - Tmin) / partitions )) for i in range(partitions) ]
listQex = [ an.getQex( 1/(listT[i] * k_Boltzmann), listEnergy, listDOS) for i in range(len(listT)) ]

plt.plot( listEnergy , log_idos, 'o-' )
# plt.plot( listEnergy , listDOS, 'o-' )
plt.plot( listT , listQex, 'o-' )




listEnergy    = [ E[i] for i in range(le)]
listDOS       = [ exp( log_idos[i] ) for i in range(le) ]
x_listT, y_listCv = an.getHeatCapacity(k_Boltzmann, Tmin, Tmax, partitions, listEnergy, listDOS)
plt.plot( x_listT, y_listCv, 'o-' )

plt.plot( x_listT, yy, 'o-' )

from math import log, exp
from scipy.interpolate import UnivariateSpline # for derivatives

logDos = [ log(listDOS[i]) for i in range(le) ]
Qex    = [ an.getQex(x_listT[i], listEnergy, listDOS) for i in range(len(x_listT))]
LogQex = [ log(Qex[i]) for i in range(len(Qex)) ]

dos_spl    = UnivariateSpline(listEnergy, listDOS, s=0, k=4)
logDos_spl = UnivariateSpline(listEnergy, logDos, s=0, k=4)
Qex_spl    = UnivariateSpline(x_listT, Qex, s=0, k=4)
LogQex_spl = UnivariateSpline(x_listT, LogQex, s=0, k=4)

plt.plot( x_listT, LogQex_spl(x_listT), 'o-' )


plt.plot( listEnergy, dos_spl(listEnergy), 'o-' )

plt.plot( listEnergy, logDos_spl(listEnergy), 'o-' )



#
