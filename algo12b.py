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



import harmonic as ha


#########################################################333
import configuration as cf
import numpy as np
from scipy.integrate import simps # Integrating using Samples
from math import log, exp
import mcFunctions5 as mc
import calc as ca

# listOfPositions = cf.getOctahedral()
# Emin = cf.getLJenergy(listOfPositions)
# forceMatrix = cf.getForceMatrix()

listOfPositions = cf.getOctahedral()
N = len(listOfPositions)
V = (15)**3 # cube of 15Angstrom as in the pymatnest tutorial
epsilon = 0.1
sigma   = 2.5
rc      = 7.50
getForceMatrix = True
Emin, eig = cf.getLJeigenvalues(listOfPositions, epsilon, sigma,rc, getForceMatrix)
Emin
eig

Emin

Emin = cf.getLJenergy(listOfPositions, sigma, epsilon)
Emin

maxDeltaE = 0.5 * sum(eig) * ( 0.01 ** 2)
maxDeltaE

deltaE = maxDeltaE #/ 100
E1 = Emin + deltaE
E2, I1, I2, errorPercent = ca.getEnext(eig, N, V, Emin, E1)
E2, I1, I2, errorPercent

# Plot
from matplotlib import pyplot as plt
intervalos = 50
x  = np.array( [ (i * deltaE / intervalos) for i in range(intervalos +8) ] )
y  = cf.harmonicDOS(x, eig, N, V) # x departs from Emin -> interpreted as dE
y1 = cf.harmonicDOS(deltaE, eig, N, V)
y2 = cf.harmonicDOS(E2-Emin, eig, N, V)


x = x + Emin
plt.plot(x, y, '-',   [Emin, E1, E2], [0,y1,y2], 'o')


x2 = np.array( [ x[i] for i in range( 20, len(x) )] )
y2 = np.array( [ log(y[i]) for i in range( 20, len(y) )] )
len(x2)
len(y2)
plt.plot(x2, y2, '-',)


################################################################################
log2 = log(2)  #math.log(2.0)


# Because python lists begin with index in zero, then `m` will begin
# in zero, instead of 1.

log_idos = [0, 0] # = [log_idos_m=0, log_idos_m=1] #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# log_idos = [0, 0] # = [log_idos_m=0, log_idos_m=1]
log_sum_idos = log(1 + 1) # = log( sum(exp(log_idos))) = log(1 + 1)

# α_m=0 = 1: impossed
# α_m=1 = 1: impossed
# α_m=2 = 1: because log_idos[0] = log_idos[1]
# lAlpha = [1, 1, 1] # [α0, α1, α2, ..., αm, ...]
lAlpha = [1, 1, 1] # [α0, α1, α2, ..., αm, ...] #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
m = 2  # python begins in 0. << the number of the next subdivision (instead of m=3 as in the paper)
#m = 3  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


################################################################################
# E_m values:
E = [Emin, E1, E2] # = [E0, E1, E2]
E

#test
E3temp, I1, I2, errorPercent = ca.getEnext(eig, N, V, Emin, E2)
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
L = 0.06 #0.01 #~amplitud of random walk
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

def draw(X):
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


#
def move(X, sigma, epsilon, nSteps, Emin, Em, Emplus1, L):
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
        Xt = copy.deepcopy(X)

        # d = Vector3(getRandomWalk(L), getRandomWalk(L), getRandomWalk(L))
        # d = unitVect(idim) * getRandomWalk(L)
        d = unitVect2(idim) * L * (random() - 0.5)
        # d = unitVect2(idim) * l



        Xt[iAtom] += d
        e  = cf.getLJenergy(Xt, sigma, epsilon)
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
            draw(X)
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
    ratioAcceptances = acceptances / nSteps
    alpha = I_EmToEmplus1 / I_EminToEm
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

def goodMove(Xo, sigma, epsilon, nSteps, Emin, Em, Emplus1, log_idos, log_sum_idos, L):
    ratioMax = 0.35
    ratioMin = 0.25
    good = False
    count = 0
    while ((not good) and count <= 10):
        count += 1
        alpha, ratioAcceptances = move(Xo, sigma, epsilon, nSteps, Emin, Em, Emplus1, L)
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
    log_idos.append( log_sum_idos + log(alpha) ) # =log_idos[m]
    DOSintegrals = [ exp(log_int) for log_int in log_idos ]
    log_sum_idos = log(sum(DOSintegrals))

    # if good:
        # log_idos.append( log_sum_idos + log(alpha) ) # =log_idos[m]

        # calculate log_sum_idos:
        # Log of the sum of integrals for each interval, m=1, m=2, m=3,...
        # DOSintegrals = [ exp(log_int) for log_int in log_idos ]
        # log_sum_idos = log(sum(DOSintegrals))

    return alpha, ratioAcceptances, log_idos, log_sum_idos, L
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
# ii
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
    alpha, ratioAcceptances, log_idos, log_sum_idos, L  =\
                goodMove(Xo, sigma, epsilon, nSteps, Emin, Em, Emplus1, log_idos, log_sum_idos, L)

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
        alpha, ratioAcceptances, log_idos, log_sum_idos, L  =\
                goodMove(Xo, sigma, epsilon, nSteps, Emin, Em, Emplus1, log_idos, log_sum_idos, L)

        # lCfgs, alpha, log_idos, log_sum_idos, L =\
        #             mc.randomMCmoves(Emin, E[m - 2], E[m - 1], Einfty,\
        #                           nSteps, e, X, log_idos, log_sum_idos, L)

        lAlpha[m] = alpha # <- overwrite since α_m was recalculated.
    #
#
X
E
################################################################################
# normalize Ω (taking advantage of logarith property when substracting):
log_idos = [log_idos[i] - log_sum_idos for i in range(len(log_idos)) ]

len(E)
len(log_idos)
len(lAlpha)

le = min(len(E), len(log_idos), len(lAlpha))

for m in range(le):
    print( E[m], log_idos[m], lAlpha[m] )
#

# l = len(E)
for i in range(le):
    if (i + 1 < le): print(E[i + 1] - E[i])

ll = le - 1
plt.plot([lAlpha[i] for i in range(ll)], 'o')

plt.plot([E[i] for i in range(ll-5)], [exp(log_idos[i]) for i in range(ll-5)], 'o')

plt.plot([E[i] for i in range(ll-5)], [log_idos[i] for i in range(ll-5)], 'o')



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
# import imp
# imp.reload(an)
import analyze as an
le = min(len(E), len(log_idos))

partitions = 1000
Tmin = 100
Tmax = 1000
k_Boltzmann   = 1


listEnergy    = [ E[i] for i in range(le)]
listDOS       = [ exp( log_idos[i] ) for i in range(le) ]
x_listT, y_listCv, yy = an.getHeatCapacity(k_Boltzmann, Tmin, Tmax, partitions, listEnergy, listDOS)
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



################################################################################
from matplotlib import pyplot as plt
from math import log, exp
# from scipy.interpolate import UnivariateSpline # for derivatives
import numpy as np, fileinput, itertools, sys
fileName = "sample.txt"
inputs = fileinput.input(files=fileName)
line_skip=0
line_end=None
interval=1
lines = itertools.islice(inputs, line_skip, line_end, interval)

lines

E = []
log_idos = []
lAlpha = []
for line in lines:
    fields = line.split()
    e, idos, a = fields
    E.append(float(e))
    log_idos.append(float(idos))
    lAlpha.append(float(a))
#
import analyze as an
le = min(len(E), len(log_idos))
len(E)
len(log_idos)

plt.plot(E, log_idos, 'o')
from scipy import interpolate
f = interpolate.interp1d(E, log_idos)
xnew = np.arange(min(E), max(E), 0.001)
ynew = f(xnew)
plt.plot(E, log_idos, 'o', xnew, ynew, '-')


from scipy.misc import derivative
xnew = np.arange(-1.2, max(E), 0.001)
# xnew = [val for val in xnew[1: len(xnew)-1]]
dy = [derivative(f, val, dx=1e-6, n=1) for val in xnew ]
plt.plot(xnew, dy)

















###
