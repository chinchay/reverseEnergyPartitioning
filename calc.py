import configuration as cf
import numpy as np
from scipy.integrate import simps

def getEnext(eig, N, V, Emin, E):
    deltaE = E - Emin
    nPart  = 1000 #50 # number of partitions

    # Getting sample: (# xMin = Emin, xMax = deltaE)
    x = np.array( [ (i * deltaE / float(nPart)) for i in range(nPart + 1) ] ) #`float()` to ensure the division will be a float (python2 resembles Fortran in divisions)
    y = cf.harmonicDOS(x, eig, N, V) # x departs from Emin -> interpreted as dE

    # integrate using Simpson method:
    I1 = simps(y, x)
    # I1

    # determine E2:
    maxN = 500 #250 # maximum number of iterations to find I2 such that
    I2 = 0
    i = 0
    x = np.array( [deltaE] )
    # x
    #
    while ( (I2 <= I1) and (i < maxN) ):
        i += 1

        # Getting sample. x[i] is the distance from Emin.
        # `i` is free to increase indef
        deltaE2 = deltaE + ( i * deltaE / nPart)
        x = np.append(x, deltaE2 )
        # x
        y = cf.harmonicDOS(x, eig, N, V) # x departs from Emin -> interpreted as dE

        # integrate using Simpson method:
        I2 = simps(y, x)

        # print(i, I1, I2)
    #
    error = 100 * ((I2 - I1) / I1)
    deltaE2 = deltaE + ( (i - 1) * deltaE / nPart)
    E2 = Emin + deltaE2
    return E2, I1, I2, error
#
