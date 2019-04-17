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

import mcFunctions as mc

# Get names and equilibrium positions of atoms:
# you will need to find first the equilibrium positions when having `V`
# For the project, MTP or DFT?? will give the equilibrium positions.
lAtomsEquilibrium, a1, a2, a3 = potential.getAtomsAtEquilibriumPositions()

# Atoms as harmonic oscillators with potential V.
# Potential energy is calculated at their equilibrium positions:
#Emin = potential.getTotalPotentialEnergy(lAtomsEquilibrium, V)
Emin = potential.get3DPotNearFirstNeighb(lAtomsEquilibrium, V, a1, a2, a3)
# from potential import VLJ2 as V2
# Emin = potential.getTotalPotentialEnergy([lAtomsEquilibrium[i].position for i in range(len(lAtomsEquilibrium))], V2)

Emin

#deltaE = 0.1 #0.1 described respect to Emin: deltaE = Energy - Emin
deltaE = abs( 0.1 * Emin / len(lAtomsEquilibrium))
deltaE
E1 = Emin + deltaE

nPart  = 50 # number of partitions
Vol = float(len(lAtomsEquilibrium)) ** (1/3)   #1.0
Vol

# Getting sample: (# xMin = Emin, xMax = deltaE)
x = np.array( [ (i * deltaE / nPart) for i in range(nPart + 1) ] )
y = functions.harmonicDOS(Vol, x) # x departs from Emin -> interpreted as dE

# integrate using Simpson method:
I1 = simps(y, x)

# determine E2:
maxN = 250 # maximum number of iterations to find I2 such that
I2 = 0
i = 0
x = np.array( [deltaE] )

#
while ( (I2 <= I1) and (i < maxN) ):
    i += 1

    # Getting sample. x[i] is the distance from Emin.
    # `i` is free to increase indef
    deltaE2 = deltaE + ( i * deltaE / nPart)
    x = np.append(x, deltaE2 )

    # x departs from Emin -> interpreted as dE
    y = functions.harmonicDOS(Vol, x)

    # integrate using Simpson method:
    I2 = simps(y, x)

    print(i, I1, I2)
#

deltaE2 = deltaE + ( (i - 1) * deltaE / nPart)
E2 = Emin + deltaE2

# print(I1, I2, Emin, E1, E2)


# Plot
from matplotlib import pyplot as plt
intervalos = 50
x = np.array( [ (i * deltaE / intervalos) for i in range(6 * intervalos) ] )
y = functions.harmonicDOS(Vol, x) # x departs from Emin -> interpreted as dE
y1 = functions.harmonicDOS(Vol, deltaE)
y2 = functions.harmonicDOS(Vol, deltaE2)

x = x + Emin
plt.plot(x, y, '-',   [Emin, E1, E2], [0,y1,y2], 'o')


##=====
import math
x2 = np.array( [ x[i] for i in range( 200, len(x) )] )
y2 = np.array( [ math.log(y[i]) for i in range( 200, len(y) )] )
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

E1 - Emin
E2 - E1

################################################################################
continuar = True
iMax   = 100 #15 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
nSteps =  20000 #1000 #5000 # 500 #50#  = nBarridos * nAtomsInConfig
alphaMin = 0.01#0.5 # the same as the paper. "At this point, the newly guessed Em is
               # deleted and then a set of random (all moves are accepted) MC
               # moves are performed... more info about how to do the last part.

# select a starting configuration:
#e, X = Emin, lAtomsEquilibrium
e, X = Emin, [lAtomsEquilibrium[i].position for i in range(len(lAtomsEquilibrium))]

# Main Reverse Energy Partitioning loop:
L = 0.01  #0.01 #0.2 #0.5 # 0.1#  ~amplitud of random walk
i = 0

import copy


while ( continuar and (i < iMax) ): # iMax allows to force stopping.
    i += 1

    # Guess the energy boundary Em
    Em = mc.getEm(lAlpha[m - 1], E[m - 1], E[m - 2])
    E.append(Em)

    # X = get
    # e = mc.getEnergyConfig(X)

    # MC sampling. Collect quantities for the next subdivision m+1.
    Xcopy = copy.deepcopy(X)
    lCfgs, alpha, log_idos, log_sum_idos, L =\
                        mc.randomMCmoves(Emin, E[m - 2], E[m - 1], E[m],\
                                      nSteps, e, Xcopy, log_idos, log_sum_idos, L, a1, a2, a3)

    lAlpha.append(alpha) # alpha = α_m

    # select a starting configuration for the next subdivision:
    assert( len(lCfgs) != 0 )
    e, X = mc.select_start_config(lCfgs, E[m - 1], E[m])

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
        lCfgs, alpha, log_idos, log_sum_idos, L =\
                    mc.randomMCmoves(Emin, E[m - 2], E[m - 1], Einfty,\
                                  nSteps, e, X, log_idos, log_sum_idos, L)

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
plt.plot(E, [exp(log_idos[i]) for i in range(len(E))], 'o')

plt.plot(E, [log_idos[i] for i in range(len(E))], 'o')

plt.plot(E, lAlpha, 'o')

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
