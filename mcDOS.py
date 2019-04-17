################################################################################
# python3
import os
print(os.getcwd())
os.chdir('/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/')
print(os.getcwd())

# import graphlab
################################################################################
# first for an harmonic oscillator:

delta =
unit =

from matplotlib.pylab import *

#=================
def belongs(x, xmin, xmax):
    return (xmin <= x) and (x <= xmax)

#=================
def myDeltaFun(x, xcenter, sigma):
    if belongs(x, xcenter - sigma, xcenter + sigma):
        return 1 / (2.0 * sigma)
    else:
        return 0

#=================
def getR(n, sigma):
    I1 = 0
    I2 = 0
    for i in range(n):
        x = random() - 0.5 # -0.5 was used to center at the origin, 0.
        y = myDeltaFun(x, 0, sigma)
        if (y > 0):
            I1 += 1
        #
        # I2 += 1

        # print(x,y)
    #
    # return I1 / I2
    return I1 / (i + 1)

sigma = 0.01
rValues = [getR(i + 1, sigma) for i in range(1000)]
plot(rValues)

#=================
def harmonicEnergy1D(x, xc, k):
    Emin = 0
    return Emin + ( 0.5 * k * ( pow(x - xc, 2) ) )

def harmonicEnergy3D(x, y, z, k):
    Emin = 0
    return Emin + ( 0.5 * k * ( pow(x,2) + pow(y,2) + pow(z,2) ) )

def getRandomCenteredZero(L): # L: amplitud
    return L * (random() - 0.5) # `0.5`:for centering around 0.

#=================
def getDOS(Energy, nSamples, sigma):
    I1 = 0
    I2 = 0
    fun = 0
    for i in range(nSamples):
        xPos = getRandomCenteredZero(4)
        yPos = getRandomCenteredZero(4)
        zPos = getRandomCenteredZero(4)

        # E = harmonicEnergy1D(xPos, 0, 1)
        E = harmonicEnergy3D(xPos, yPos, zPos, 1)

        y = myDeltaFun(E, Energy, sigma)
        if (y > 0):
            # I1 += 1
            fun += y
        #
        # I2 += 1

        # print(x,y)
    #
    # return I1 / I2
    # return I1 / nSamples
    return fun / nSamples

#Energy = 0.05
#sigma = 0.01

Energy = 0.01
sigma = 0.0001

# just to get an improved value:
dos = average( [getDOS(Energy, 1000, sigma) for i in range(10)]  )
dos

n = 30
#sigma = 0.001
sigma = 0.01
Eval = []
dosVal = []
for i in range(n):
    # E = 0.05 * i / n
    E = 10 * i / n
    dos = average( [getDOS(E, 1000, sigma) for i in range(10)]  )
    Eval.append(E)
    dosVal.append(dos)
#
plot(Eval, dosVal, 'o')







myList = [i for i in range(10, 1000, 50)]
sigmas = [1/(i + 1) for i in myList]

# dosAtSigmas = [ getDOS(Energy, 1000, s) for s in sigmas ]
# plot(myList, dosAtSigmas, 'o')

avgDOSatSigmas = [ average( [getDOS(Energy, 1000, s) for i in range(10)]  ) for s in sigmas ]
plot(myList, avgDOSatSigmas, 'o')

###############################################################################
plot(rValues)





#
