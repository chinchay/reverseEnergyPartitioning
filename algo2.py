################################################################################
# python3
import os
print(os.getcwd())
os.chdir('/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/')
print(os.getcwd())

# import graphlab
################################################################################

import functions

import potential
from potential import VLJ as V

import numpy as np
from scipy.integrate import simps # Integrating using Samples

# import mcFunctions as mc
import mcFunctions2 as mc2
import harmonic as ha

import crystal as cr
import numpy as np



# Get names and equilibrium positions of atoms:
# you will need to find first the equilibrium positions when having `V`
# For the project, MTP or DFT?? will give the equilibrium positions.
# lAtomsEquilibrium, a1, a2, a3 = potential.getAtomsAtEquilibriumPositions()

# Xeq, a1, a2, a3, Emin, forceMatrix = ha.getConfigEquilibrium()
# structEq = cr.getStruct() # to be used from cfg of mtp
structEq = cr.getFccEquilibrium()
radio    = 1.1 * structEq[0].distance(structEq[1]) # for first nearest neighbors
Emin     = cr.getEnergy(structEq, radio, [])

indxNeigsFromEachSite = cr.getIndxNeigsFromEachSite(structEq, radio)
dRneighborsFromEachSite = cr.getDRneighborsFromEachSite(structEq, radio)
forceMatrix = cr.getForceMatrix(structEq, indxNeigsFromEachSite, dRneighborsFromEachSite, radio)

Emin
# forceMatrix

# structEq.num_sites
#
# structEq[0].coords
# structEq[0].y
# # structEq[0].coords + structEq[0].coords
# # structEq[0].coords = structEq[0].coords + np.array([10,20,30])
#
# structEq[1].specie
# structEq[0] = structEq[1].specie, structEq[0].coords + np.array([10,20,30])
#
# structEq[0]
#
# structEq[1].frac_coords
# ss = structEq[1]
# aa = ss.specie
# aa
# structEq[2]
# structEq.cart_coords
# i
# lAtomsEquilibrium, a1, a2, a3, Emin = ha.getAtomsAtEquilibriumPositions()



# Atoms as harmonic oscillators with potential V.
# Potential energy is calculated at their equilibrium positions:
    #Emin = potential.getTotalPotentialEnergy(lAtomsEquilibrium, V)
# Emin = potential.get3DPotNearFirstNeighb(lAtomsEquilibrium, V, a1, a2, a3)
    # from potential import VLJ2 as V2
    # Emin = potential.getTotalPotentialEnergy([lAtomsEquilibrium[i].position for i in range(len(lAtomsEquilibrium))], V2)


#deltaE = 0.1 #0.1 described respect to Emin: deltaE = Energy - Emin
deltaE = 1 #0.1  #1 == f*d
# deltaE = abs( 2.5 * Emin / len(lAtomsEquilibrium))
deltaE
E1 = Emin + deltaE


nPart  = 1000 #50 # number of partitions
# Vol = float(len(lAtomsEquilibrium)) ** (1/3)   #1.0
# Vol

# Getting sample: (# xMin = Emin, xMax = deltaE)
x = np.array( [ (i * deltaE / nPart) for i in range(nPart + 1) ] )
# x

# ha.harmonicDOS(-8)

# y = functions.harmonicDOS(Vol, x) # x departs from Emin -> interpreted as dE
y = ha.harmonicDOS(x, forceMatrix) # x departs from Emin -> interpreted as dE
# y


# integrate using Simpson method:
I1 = simps(y, x)
I1

# determine E2:
maxN = 500 #250 # maximum number of iterations to find I2 such that
I2 = 0
i = 0
x = np.array( [deltaE] )
x
#
while ( (I2 <= I1) and (i < maxN) ):
    i += 1

    # Getting sample. x[i] is the distance from Emin.
    # `i` is free to increase indef
    deltaE2 = deltaE + ( i * deltaE / nPart)
    x = np.append(x, deltaE2 )
    x

    # x departs from Emin -> interpreted as dE
    # y = functions.harmonicDOS(Vol, x)
    y = ha.harmonicDOS(x, forceMatrix)

    # integrate using Simpson method:
    I2 = simps(y, x)

    print(i, I1, I2)
#

deltaE2 = deltaE + ( (i - 1) * deltaE / nPart)
deltaE2 - deltaE


E2 = Emin + deltaE2
deltaE2
# print(I1, I2, Emin, E1, E2)


# Plot
from matplotlib import pyplot as plt
intervalos = 50
x = np.array( [ (i * deltaE / intervalos) for i in range(intervalos +3) ] )
# y = functions.harmonicDOS(Vol, x) # x departs from Emin -> interpreted as dE
# y1 = functions.harmonicDOS(Vol, deltaE)
# y2 = functions.harmonicDOS(Vol, deltaE2)
y  = ha.harmonicDOS(x, forceMatrix) # x departs from Emin -> interpreted as dE
y1 = ha.harmonicDOS(deltaE, forceMatrix)
y2 = ha.harmonicDOS(deltaE2, forceMatrix)


x = x + Emin
plt.plot(x, y, '-',   [Emin, E1, E2], [0,y1,y2], 'o')


##=====
# import math
# x2 = np.array( [ x[i] for i in range( 200, len(x) )] )
# y2 = np.array( [ math.log(y[i]) for i in range( 200, len(y) )] )
# len(x2)
# len(y2)
# plt.plot(x2, y2, '-',)
import math
x2 = np.array( [ x[i] for i in range( 20, len(x) )] )
y2 = np.array( [ math.log(y[i]) for i in range( 20, len(y) )] )
len(x2)
len(y2)
plt.plot(x2, y2, '-',)

# Ln( DOS) is approx. linear respect to Energy!
# this is true for D > 2 (for electrons nao interagentes temos D=1 e a
# densidade de estados é proporcional a E^(-1/2), e a DOS = segunda Derivada de
# N con respecto a E da negativo. Porem, se tivessemos mais graus de liberdade,
# que é este caso de atoms interagindo em um cristal? temo segunda derivada
# positiva)

# Still, analyze carefully the interaction, the matrix force, are we dealing
# with a crystal? or with a finite number of atoms?

#=================
# already defined in functions.py:
# def belongs(x, xmin, xmax):
#     return (xmin <= x) and (x <= xmax)

from math import log, exp
log2 = log(2)  #math.log(2.0)


# Because python lists begin with index in zero, then `m` will begin
# in zero, instead of 1.

log_idos = [0, 0, 0] # = [log_idos_m=0, log_idos_m=1] #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# log_idos = [0, 0] # = [log_idos_m=0, log_idos_m=1]
log_sum_idos = log(1 + 1) # = log( sum(exp(log_idos))) = log(1 + 1)

# α_m=0 = 1: impossed
# α_m=1 = 1: because log_idos[0] = log_idos[1]
# α_m=2 = ? to find for m >=2 (in the paper `m` will be 3)
# lAlpha = [1, 1] # [α0, α1, α2, ..., αm, ...]
lAlpha = [1, 1, 1] # [α0, α1, α2, ..., αm, ...] #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# m = 2 # the number of the next subdivision (instead of m=3 as in the paper)
m = 3  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# E_m values:
E = [Emin, E1, E2] # = [E0, E1, E2]

E


[Emin - Emin, E1 - Emin, E2 - Emin, mc2.getEm(lAlpha[m - 1], E[m - 1], E[m - 2]) - Emin]

[E1-Emin, E2-E1, mc2.getEm(1, E[3 - 1], E[3 - 2]) -E2]
[E1-Emin, E2-E1, E2+(E2-E1)/(E1-Emin) -E2]

mc2.getWeigth(E2 + 0.0001, E[m - 1], E[m - 2])



################################################################################
continuar = True
iMax   = 100 #30 #10 #100 #15 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
nSteps =  20000 #5000 #20000 #1000 #5000 # 500 #50#  = nBarridos * nAtomsInConfig
alphaMin = 0.01#0.5 # the same as the paper. "At this point, the newly guessed Em is
               # deleted and then a set of random (all moves are accepted) MC
               # moves are performed... more info about how to do the last part.

# select a starting configuration:
#e, X = Emin, lAtomsEquilibrium
# e, X = Emin, [lAtomsEquilibrium[i].position for i in range(len(lAtomsEquilibrium))]

import copy
e = Emin
# X = copy.deepcopy(Xeq)
structX = copy.deepcopy(structEq)
# structX
# Main Reverse Energy Partitioning loop:
#L = 0.005 #0.05 #0.1 #0.01  #0.01 #0.2 #0.5 # 0.1#  ~amplitud of random walk
L = 0.01 #0.01

L0 = L
i = 0
nSteps0 = nSteps

# aa = [10,20, 20,30,40]
# aa.index(20)
# aa = [6,5,4,4,2,1,10,1,2,48]
# set(aa)
# sorted(set(aa))



import mcFunctions2 as mc2
Em = mc2.getEm(lAlpha[m - 1], E[m - 1], E[m - 2])
Em
E.append(Em)
structXcopy = copy.deepcopy(structX)
mc2.randomMCmoves(Emin, E[m - 2], E[m - 1], E[m],\
              nSteps, e, structXcopy, log_idos, log_sum_idos, L, structEq, forceMatrix, radio, dRneighborsFromEachSite)


iii

while ( continuar and (i < iMax) ): # iMax allows to force stopping.
    i += 1

    # # Guess the energy boundary Em
    # if i == 1:
    #     Em = E2 + (E2-E1)/(E1-Emin)
    # else:
    #     Em = mc.getEm(lAlpha[m - 1], E[m - 1], E[m - 2])

    Em = mc2.getEm(lAlpha[m - 1], E[m - 1], E[m - 2])
    Em
    E.append(Em)


    [Emin-Emin, E[m - 2]-Emin, E[m - 1]-Emin, E[m]-Emin]

    # X = get
    # e = mc.getEnergyConfig(X)

    # MC sampling. Collect quantities for the next subdivision m+1.
    structXcopy = copy.deepcopy(structX)
    # lCfgs, alpha, log_idos, log_sum_idos, L =\
    #                     mc.randomMCmoves(Emin, E[m - 2], E[m - 1], E[m],\
    #                                   nSteps, e, Xcopy, log_idos, log_sum_idos, L, a1, a2, a3, Xeq, forceMatrix)
    lCfgs, alpha, log_idos, log_sum_idos, L =\
                        mc2.randomMCmoves(Emin, E[m - 2], E[m - 1], E[m],\
                                      nSteps, e, structXcopy, log_idos, log_sum_idos, L, a1, a2, a3, structEq, forceMatrix, radio)

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



    L = L0
    nSteps = nSteps0
    lAlpha.append(alpha) # alpha = α_m

    # select a starting configuration for the next subdivision:
    assert( len(lCfgs) != 0 )
    nthMin = 1
    e, X = mc.select_start_config(lCfgs, E[m - 1], E[m], nthMin)

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
        # lCfgs, alpha, log_idos, log_sum_idos, L =\
        #             mc.randomMCmoves(Emin, E[m - 2], E[m - 1], Einfty,\
        #                           nSteps, e, X, log_idos, log_sum_idos, L, a1, a2, a3, Xeq, forceMatrix)
        lCfgs, alpha, log_idos, log_sum_idos, L =\
                    mc2.randomMCmoves2(Emin, E[m - 2], E[m - 1], Einfty,\
                                  nSteps, e, structXcopy, log_idos, log_sum_idos, L, a1, a2, a3, structEq, forceMatrix)

        lAlpha[m] = alpha # <- overwrite since α_m was recalculated.
    #
#
E

# normalize Ω (taking advantage of logarith property when substracting):
log_idos = [log_idos[i] - log_sum_idos for i in range(len(log_idos)) ]

for m in range(len(E)):
    print( E[m], log_idos[m], lAlpha[m] )
#

l = len(E)
for i in range(l):
    if (i + 1 < l): print(E[i + 1] - E[i])

len(log_idos)
len(E)
nn = len(E) #- 2
plt.plot([E[i] for i in range(nn)], [exp(log_idos[i]) for i in range(nn)], 'o')

plt.plot([E[i] for i in range(nn)], [log_idos[i] for i in range(nn)], 'o')

plt.plot([E[i] for i in range(nn)],  [lAlpha[i] for i in range(nn)], 'o')

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
