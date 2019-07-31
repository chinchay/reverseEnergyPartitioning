from scipy.integrate import simps # for simps integration
from math import exp, log
# from scipy.interpolate import UnivariateSpline # for derivatives
import numpy as np
from scipy import interpolate
from scipy.misc import derivative

def intDOSinEnergy(listEnergy, listDOS):
    assert(len(listEnergy) == len(listDOS))

    # Getting sample: (# xMin = Emin, xMax = deltaE)
    x = np.array( listEnergy )
    y = np.array( listDOS )

    # integrate using Simpson method:
    return simps(y, x)
#


def getQex(beta, listEnergy, listDOS):
    assert(len(listEnergy) == len(listDOS))

    Qex = 0
    for i in range(len(listDOS)):
        Qex += listDOS[i] * exp(-listEnergy[i] * beta)
    #
    Qex = Qex / intDOSinEnergy(listEnergy, listDOS)
    return Qex
#

def getHelmholtzFreeEnergy(beta):
    return -log( getQex(beta) ) / beta
#

def getHeatCapacity(k_Boltzmann, Tmin, Tmax, partitions,  listEnergy, listDOS):
    # https://stackoverflow.com/questions/40226357/second-derivative-in-python-scipy-numpy-pandas
    # get second derivative of the interpolated function:

    betaMin = 1.0 / (k_Boltzmann * Tmax)
    betaMax = 1.0 / (k_Boltzmann * Tmin)
    n = partitions

    x_listBeta = [ betaMin + (i * ((betaMax - betaMin) / n) ) for i in range(n) ]
    y_LnQex    = [ log( getQex( x_listBeta[i], listEnergy, listDOS ) ) for i in range(n) ]

    f          = interpolate.interp1d(x_listBeta, y_LnQex)
    # D2_f       = [ derivative(f, val, dx=1e-6, n=2) for val in x_listBeta ]
    # y_spl_2d   = y_spl.derivative(n=2)

    ###
    x_listT    = [ Tmin + (i * ( (Tmax - Tmin) / n )) for i in range(1,n-1) ]
    listBeta   = [ 1 / (x_listT[i] * k_Boltzmann) for i in range(len(x_listT)) ]
    D2_f       = [ derivative(f, val, dx=1e-6, n=2) for val in listBeta ]
    y_listCv   = [ k_Boltzmann * (listBeta[i] ** 2) * D2_f[i] for i in range(len(x_listT)) ]

    ###
    # listQex    = [ getQex( listBeta[i], listEnergy, listDOS ) for i in range(n) ]
    # y_Qex = UnivariateSpline(x_listBeta, y_LnQex, s=0, k=4)

    # plt.plot( x_listT, y_listCv )
    # y_Qex = UnivariateSpline(x_listT, [], s=0, k=4)
    return x_listT, y_listCv#, y_spl(listBeta),
#
#
