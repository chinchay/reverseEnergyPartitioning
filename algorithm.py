

# Set E1, E2 using the harmonic approximation

import functions

import potential
from potential import VLJ as V

import numpy as np
from scipy.integrate import simps # Integrating using Samples

# Get names and equilibrium positions of atoms:
# you will need to find first the equilibrium positions when having `V`
# For the project, MTP or DFT?? will give the equilibrium positions.
lAtomsEquilibrium = potential.getAtomsAtEquilibriumPositions()

# Atoms as harmonic oscillators with potential V.
# Potential energy is calculated at their equilibrium positions:
Emin = potential.getTotalPotentialEnergy(lAtomsEquilibrium, V)

deltaE = 0.1 # described respect to Emin: deltaE = Energy - Emin
E1 = Emin + deltaE

nPart  = 50 # number of partitions
Vol = 1.0

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



log_idos = [0]
log_idos.append(0)

# α0 = 1: impossed
# α1 = 1: impossed
# α(m=2): because log_idos[0] = log_idos[1]
# αm = ? to find for m >=3
α = [1, 1, 1] # [α0, α1, α2, ..., αm, ...]
m = 3 # the number of the next subdivision

# Em values:
E = [Emin, E1, E2]

log2 = math.log(2.0)

# Main Reverse Energy Partitioning loop:

# Guess the energy boundary Em
minimo = min( [ α[m - 1], 1 ] )
rLog = math.log( 1 +  1 / minimo ) / log2
Em = E[m - 1] + (  (E[m - 1] - E[m - 1]) * rLog  )
