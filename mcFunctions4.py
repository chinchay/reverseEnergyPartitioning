#call mc_move():
    #get Xold -> Xnew, small enough such that the `E` resultant does not go outside of [Emin,Em]
    #accept the movement Xold->Xnew with probability P(o->n)=min(1, w(En)/w(Eo))
    #get the new energy of Xnew: E

from math import log, exp, cos, sin, pi
# log2 = log(2)
import os


from functions import belongs
from matplotlib.pylab import *  # for random()

def getEm(α_m, E_m, E_mMinus1):
    # minimo = min( [ α[m - 1], 1 ] )
    # rLog   = log( 1 +  (1 / minimo) ) / log2
    # Em     = E[m - 1] + (  (E[m - 1] - E[m - 2]) * rLog  )
    minimo = min([ 1, α_m ])
    # log3 = log2 + 1.0
    rLog   = log( 1 +  (1 / minimo) ) / log(2) #<<<<<<<<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<<<<<<<
    return E_m + (  (E_m - E_mMinus1) * rLog  )
#

def getRandomWalk(L): # L: amplitud
    return L * (random() - 0.5) # `0.5`:for centering around 0.

def myDiscreteRandom():
    if random() > 0.5 :
        return 1
    else:
        return 0
#

def myDiscreteRandom3():
    if belongs(random(),0, 0.3333):
        return 1
    elif belongs(random(),0.3333, 0.6666):
        return 2
    else:
        return 3
#

# def myDiscreteRandom6():
    # n = 6
    # r = random()
    # if belongs(r, 0, 1.0/n):
    #     return 1
    # elif belongs(r, 1.0/n, 2/n):
    #     return 2
    # elif belongs(r, 2.0/n, 3.0/n):
    #     return 3
    # elif belongs(r, 3.0/n, 4.0/n):
    #     return 4
    # elif belongs(r, 4.0/n, 5.0/n):
    #     return 5
    # else:
    #     return 6
    # #
    #
#


def getRandomWalk3D(L):
    # L = 0.1 # think a better factor <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!!
    # L is related to ~diffusion coefficient??? L must be small enough

    # ii = myDiscreteRandom6()


    ii = int( 6 * random() )

    l = L
    # l = getRandomWalk(L)
    if ii == 0:
        delta = [0, 0, -l]
    elif ii == 1:
        delta = [0, 0,  l]
    elif ii == 2:
        delta = [0, -l, 0]
    elif ii == 3:
        delta = [0,  l, 0]
    elif ii == 4:
        delta = [-l, 0, 0]
    else:
        delta = [ l, 0, 0]
    #
    return delta





    # dim = int( random() * 3 )
    # dim = myDiscreteRandom3()
    # if dim == 1:
    #     return [getRandomWalk(L), 0, 0] # = [dx,0,0]
    # elif dim == 2:
    #     return [0, getRandomWalk(L), 0] # = [0, dy, 0]
    # elif dim == 3:
    #     return [0, 0, getRandomWalk(L)] # = [0,0,dz]

    # r     = random() * L
    # theta = random() * pi
    # phi   = random() * 2 * pi
    # x = r * cos(theta) * cos(phi)
    # y = r * cos(theta) * sin(phi)
    # z = r * sin(theta)
    # return [x, y, z]

    # return [L, getRandomWalk(L), getRandomWalk(L)] # = [dx,dy,dz]
    # return [getRandomWalk(L), getRandomWalk(L), getRandomWalk(L)] # = [dx,dy,dz]
    # return [getRandomWalk(L), getRandomWalk(L), 0] # = [dx,dy,dz] <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#



def addVectors(a, b): # `a` and `b` must have the same lenght!
    from vectormath import Vector3 as Vector3
    return Vector3([a[i] + b[i] for i in range(len(a))])



def getNewConfig(i, X, L): # L ~amplitud of random walk
    # X_old is a vector with [r1,r2,r3,...,rNatoms] with the positions of every
    # atom of the configuration. ri = [x,y,z]
    iAtom = i % len(X)
    # iAtom = int( random() * (len(X) - 1) ) # << let's avoid to move the last atom to be taken as reference
    # iAtom = int( random() * len(X) )

    X[iAtom] = addVectors( X[iAtom], getRandomWalk3D(L) )

    # for j in range(len(X)-1): # << let's avoid to move the last atom to be taken as reference
        # X[j] = addVectors( X[j], getRandomWalk3D(L) )
        # if myDiscreteRandom() == 1:
            # X[j] = addVectors( X[j], getRandomWalk3D(L) )
    #


    return X
#

################################################################################
def getEnergyConfig(X):
    import potential
    from potential import VLJ as V
    import harmonic as ha
    # E = 0

    import configuration as cf
    epsilon = 0.1
    sigma   = 2.5
    rc      = 7.50
    getForceMatrix = False
    E, eig = cf.getLJeigenvalues(X, epsilon, sigma, rc, getForceMatrix)
    return E
#
################################################################################


def getWeigth(E, E_m, E_m_minus_1):
    # see the correct formula in https://pubs.acs.org/doi/ipdf/10.1021/ct3007056
    if (E_m <= E):
        f = -3 * log(2)
        g = (E - E_m) / (E_m - E_m_minus_1)
        return exp(f * g)
    else: # (E < E_m_minus_1)
        return 1.0
#

def getWeigth2(E, E_mplus1):
    # see the correct formula in https://pubs.acs.org/doi/ipdf/10.1021/ct3007056
    if (E <= E_mplus1):
        return 1
    else:
        return 0
#

def getProbTransition(E_new, E_old, E_m, E_m_minus_1):
    w_new = getWeigth(E_new, E_m, E_m_minus_1)
    w_old = getWeigth(E_old, E_m, E_m_minus_1)
    return min(1, w_new / w_old)

def getProbTransition2(E_new, E_old, E_mplus1):
    w_new = getWeigth2(E_new, E_mplus1)
    w_old = getWeigth2(E_old, E_mplus1)
    assert(w_old != 0) # it should always be less than Emplus1 since it was never accepted to a high energy.
    return min(1, w_new / w_old)

# Static variable, for counting:
def counter():
    if 'cnt' not in counter.__dict__:
        counter.cnt = 0
    counter.cnt += 1
    return counter.cnt

def counterMin():
    if 'cnt' not in counterMin.__dict__:
        counterMin.cnt = 0
    counterMin.cnt += 1
    return counterMin.cnt

# def moveHarmonic(i, X, Xeq, forceMatrix, L):
#     import harmonic as ha
#     X = getNewConfig(i, X, L)
#     rHyper, deltaEharmonic = ha.getHarmonicEnergy(X, Xeq, forceMatrix)
#     return rHyper, deltaEharmonic

def mc_move(i, E_old, X_old, E_mplus1, E_m, E_m_minus_1, L): #Xeq,forceMatrix
    import copy
    import harmonic as ha

    Xtemp  = copy.deepcopy(X_old)
    X_temp = getNewConfig(i, Xtemp, L)
    X_new  = copy.deepcopy(X_temp)
    E_new  = getEnergyConfig(X_new)

    # P_old2new = getProbTransition(E_new, E_old, E_m, E_m_minus_1)
    # P_old2new = getProbTransition2(E_new, E_old, E_mplus1)


    hasMoved = False

    # if ( P_old2new >= random() ):
    if E_new <= E_mplus1: #P_old2new = 1
        hasMoved = True

        ########################################################################################
        import os
        c = counter()
        # if (c % 5 == 0):
        # if (c % 100 == 0):
        if (c % 1000 == 0):
            nAtoms = len(X_old)
            out = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/moves." + str(c) + ".xyz"
            pos = str(nAtoms) + "\n\n"
            for i in range(len(X_new)):
                pos += "H   " + str(X_new[i][0]) + "   " + str(X_new[i][1]) + "   " + str(X_new[i][2]) + "\n"
            #
            pos += "\n"

            f = open(out, "w")
            f.write(pos)
            f.close()

            #######################################################
            # rHyper, deltaEharmonic = ha.getHarmonicEnergy(X_new, Xeq, forceMatrix)
            # out = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/rHyperEHarmonic" + str(c) + ".txt"
            # s = str(rHyper) + " , " + str(deltaEharmonic) + " , " +  str(E_new) + "\n"
            # # s = str(rHyper) + " , " + str(deltaEharmonic) + "\n"
            # f = open(out, "w")
            # f.write(s)
            # f.close()


        return E_new, X_new, hasMoved
    #
    # X = copy.deepcopy(X_old)
    return E_old, X_old, hasMoved

    # e, X, hasMoved = mc_move(i, e, X, Xeq, forceMatrix, E_mMinus1, E_mMinus2, L, a1, a2, a3)

    # c = counter()
    # rHyper, deltaEharmonic = moveHarmonic(i, X_old, Xeq, forceMatrix, L)
    # out = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/rHyperEHarmonic" + str(c) + ".txt"
    # # s = str(rHyper) + " , " + str(deltaEharmonic) + " , " +  str(E_new) + "\n"
    # s = str(rHyper) + " , " + str(deltaEharmonic) + "\n"
    # f = open(out, "w")
    # f.write(s)
    # f.close()


    # X_new = copy.deepcopy( X_old )
    # X_new = getNewConfig(i, X_new, L) # L ~amplitud of random walk
    # E_new, energyPerParticle = getEnergyConfig(X_new, a1, a2, a3)
    # P_old2new = getProbTransition(E_new, E_old, E_m_minus_1, E_m_minus_2)
    # hasMoved = False

#
#
#
#     if P_old2new >= random():
#         X_old = copy.deepcopy(X_new)
#         E_old = E_new
#         hasMoved = True
#         ########################################################################################
#         import os
#         c = counter()
#         if (c % 5 == 0):
#         # if (c % 1000 == 0):
#             nAtoms = len(X_old)
#             out = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/moves." + str(c) + ".xyz"
#             pos = str(nAtoms) + "\n\n"
#             for i in range(len(X_new)):
#                 pos += "H   " + str(X_new[i][0]) + "   " + str(X_new[i][1]) + "   " + str(X_new[i][2]) + "\n"
#             #
#             pos += "\n"
#
#             f = open(out, "w")
#             f.write(pos)
#             f.close()
#
#             ########################################################
#             rHyper, deltaEharmonic = ha.getHarmonicEnergy(X_old, Xeq, forceMatrix)
#             out = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/rHyperEHarmonic" + str(c) + ".txt"
#             # s = str(rHyper) + " , " + str(deltaEharmonic) + " , " +  str(E_new) + "\n"
#             s = str(rHyper) + " , " + str(deltaEharmonic) + "\n"
#             f = open(out, "w")
#             f.write(s)
#             f.close()
#
#         #
#         return E_old, X_old, hasMoved
#
#         #     # also, save the distances respect to the bottom of Lennard-Jones potential:
#         #     eqFile = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/equilibriumPositions.txt"
#         #     pos = str(nAtoms) + "\n\n"
#         #     with open(eqFile , "r") as ifile:
#         #         i = 0
#         #         for line in ifile:
#         #             Xequil = [ float( line.split()[j] ) for j in range(3) ]
#         #             dx = Xequil[0] - X_new[i][0]
#         #             dy = Xequil[1] - X_new[i][1]
#         #             dz = Xequil[2] - X_new[i][2]
#         #             # print(dx, dy, dz)
#         #             deltaR = ( dx**2 + dy**2 + dz**2 ) ** 0.5
#         #             pos += "H   " + str(deltaR) + "   " + str(energyPerParticle[i]) + "   " + str(0) + "\n"
#         #             i += 1
#         #         #
#         #         pos += "\n"
#         #     #
#         #     # feqFile.close()
#         #     outDeltaR = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/deltaRmoves." + str(c) + ".xyz"
#         #     foutDeltaR = open(outDeltaR,'w')
#         #     foutDeltaR.write(pos)
#         #     foutDeltaR.close()
#         # ########################################################################################
#
#     #
#     # print(E_old, E_new, P_old2new,  hasMoved, E_m_minus_1, E_m_minus_2, L)
#     return E_old, X_old, hasMoved
#     # return E_old, X_old, hasMoved
# #


#

################################################################################
# MC sampling with weight w(E). Collect quantities for the next subdivision m+1.
################################################################################
def randomMCmoves(Eminimo, E_mMinus1, E_m, E_mplus1,\
                  nSteps, e, X, log_idos, log_sum_idos, L):
    import copy
    import os
    # import harmonic as ha

    Eo  = getEnergyConfig(X)
    Xo = copy.deepcopy(X)

    lCfgs = [] # to save cfgs for the next subdivision.
    # where: lCfgs_m = [cfg_1, cfg_2, ..., cfg_i, ...]
    # where: cfg_i = [energy_i, X_i ]
    # where: cfg_i = [energy_i, X_i ]
    # where: X_i = [r1, r2, r3, r4, ..., rNatoms]
    # where: r_n = [x,y,z]

    #===========================================================================
    # paper: "The maximum displacement of the translational moves is
    # automatically changed after each partitioning process to adjust the
    # acceptance toward the range 25%–35%. If the acceptance rate during a
    # partitioning process falls below 10% or above 70%, the partitioning
    # process is repeated."
    # L ~amplitud of random walk of a ...
    #             ...particle (dx,dy,dz) = (L*random(), L*random(), L*random())
    ratioOfCfgsToSave = 0.5
    ratioAcceptances  = 0
    # L      = 0.1
    repMax = 10 #10 # maximum number of repetitions to get the appropriate ratio of
               # acceptances at a certain L value.
    repeat = 0

    ratioMin = 0.25
    ratioMax = 0.35

    while( (not belongs(ratioAcceptances, ratioMin, ratioMax) and\
            (repeat < repMax) ) or\
            continuar):

        # print(condition1, condition2, continuar, finalCondition)
        continuar = False
        # Reset energy histogram: ehist = [I1, I2].
        # I1 = c*∫_Em-1^Em Ω(E)dE;  I2 = c*∫_Emin^Em-1 Ω(E)dE, `c` is some constant.


        if (repeat > 0):
            # if I2 == 0:
            #     repeat = 0
            #     # L = 0.99 * L
            if ratioAcceptances > ratioMax:
                L = 1.1 * L
                print("...... increasing L = ", float("{0:.3f}".format(L)), repeat, ratioAcceptances)
                X = copy.deepcopy(Xo)
                return 0
            elif ratioAcceptances < ratioMin:
                L = 0.9 * L
                print("...... decreasing L = ", float("{0:.3f}".format(L)), repeat, ratioAcceptances)
                X = copy.deepcopy(Xo)
                return 0
        #

        [I_EminToEm, I_EmToEmplus1, I_Emplus1ToInfty] = [0, 0, 0]
        repeat += 1
        acceptances = 0

        # L = L * 0.8 # Decrease the amplitud of the random walk to get the
                    # desired ratioAcceptances range values [25%, 35%].
                    # You can improve with a more sophisticated algorithm here.

        # Collect ehist for [Emin,E[m-1]] and [E[m-1],Em]:
        for i in range(nSteps):
            # get a configuration e,X after randomly move X:




            Eold = getEnergyConfig(X)
            # assert( Eo == getEnergyConfig(X) )
            hasMoved = False
            Xtemp = copy.deepcopy(X)
            Xtemp = getNewConfig(i, Xtemp, L)
            Enew  = getEnergyConfig(Xtemp)
            # P_old2new = getProbTransition(E_new, E_old, E_m, E_m_minus_1)
            # P_old2new = getProbTransition(Enew, Eold, E_mplus1, E_m)

            os.system("echo "+ str(Enew) + " >> out2.txt")

            # if (Enew <= E_mplus1 + 0.001):
            if (Enew <= E_mplus1 + 0.0006):
            # if ( P_old2new >= random() ):
                os.system("echo "+ str(Enew) + " >> out.txt")
                hasMoved = True
                e = Enew
                X = copy.deepcopy(Xtemp)
            #
            # hasMoved = True
            # e = Enew
            # X = copy.deepcopy(Xtemp)



            # e, X, hasMoved = mc_move(i, e, X, E_mplus1, E_m, E_mMinus1, L)

            # # Xo = copy.deepcopy(X)
            # # assert(id(Xo) != id(X))
            # # e, Xtemporal, hasMoved = mc_move(i, e, Xo, Xeq, forceMatrix, E_mMinus1, E_mMinus2, L, a1, a2, a3)
            #
            # # assert(id(Xo) != id(Xtemporal))
            # # if (hasMoved == False):
            # #     for k in range(len(Xo)):
            # #         for l in range(3) :
            # #             assert(Xo[k][l] == Xtemporal[k][l])
            # #         #
            # # X = copy.deepcopy(Xtemporal)
            #
            #
            # c = counter()
            # rHyper, deltaEharmonic = ha.getHarmonicEnergy(X, Xeq, forceMatrix)
            # out = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/rHyperEHarmonic" + str(c) + ".txt"
            # # s = str(rHyper) + " , " + str(deltaEharmonic) + " , " +  str(E_new) + "\n"
            # s = str(rHyper) + " , " + str(deltaEharmonic) + "\n"
            # f = open(out, "w")
            # f.write(s)
            # f.close()



            # it has no sense to add the same configuration (not moved) to Ω:
            if hasMoved:
                acceptances += 1 # accepted move.
                if belongs(e, Eminimo, E_m): # (Emin <= e <= E[m - 1]): #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    I_EminToEm += 1  # un gol más para ∫_Emin^Em-1 Ω(E)dE
                elif belongs(e, E_m, E_mplus1): # (E[m - 1] <= e <= Em):
                    I_EmToEmplus1 += 1 # un gol más para ∫_Em-1^Em Ω(E)dE
                elif E_mplus1 < e:
                    I_Emplus1ToInfty += 1
                #

                # randomly save some configs to use for the next subdivision:
                if random() > ratioOfCfgsToSave: lCfgs.append([e, X])
            #
            if (e < 1.1 * Eminimo):
                print(e, Eminimo)
                ########################################################################################
                import os
                c = counter()
                nAtoms = len(X)
                out = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/moves." + str(c) + ".xyz"
                pos = str(nAtoms) + "\n\n"
                for i in range(len(X)):
                    pos += "H   " + str(X[i][0]) + "  " + str(X[i][1]) + str(X[i][2]) + "\n"
                #
                pos += "\n"

                f = open(out, "w")
                f.write(pos)
                f.close()
                ########################################################################################
                assert (1 == 0)
            # print(I1, I2, e, hasMoved, belongs(e, E_mMinus1 - 10, E_mMinus1), belongs(e, Eminimo, E_mMinus1), belongs(e,E_mMinus1, E_m), L)
        #
        ratioAcceptances = acceptances / nSteps
        # print(ratioAcceptances, I_EminToEm, I_EmToEmplus1, e, E_m, E_mplus1, belongs(e, Eminimo, E_m), belongs(e,E_m, E_mplus1), L, repeat)
        # print(ratioAcceptances, I_EminToEm, I_EmToEmplus1, I_Emplus1ToInfty, L, repeat)
        f  = (I_EminToEm + I_EmToEmplus1) / (I_EminToEm + I_EmToEmplus1 + I_Emplus1ToInfty)
        ff = I_EmToEmplus1 / (I_EmToEmplus1 + I_Emplus1ToInfty)
        print( float("{0:.3f}".format(ratioAcceptances)), I_EminToEm, I_EmToEmplus1, I_Emplus1ToInfty, "...", float("{0:.3f}".format(f)), float("{0:.3f}".format(ff)), L, repeat, Eo)

        factor =  (I_EmToEmplus1 + 1) / (I_EminToEm + 1) # '+1' to avoid zeros
        assert( len(lCfgs) != 0 )
        # if ( (I_EmToEmplus1 < 100) or (I_EminToEm < 100)):
        #     # choose X with lower energy
        #     print("A: choosing config with min energy ...")
        #     nthMin = 1
        #     e, X = select_start_config(lCfgs, E_mMinus1, E_m, nthMin)
        #     continuar = True
        if factor < 1/100.0:
            # # choose X with higher energy
            print("B: choosing config with high energy ...")
            # nthMin = 2
            # e, X = select_start_config(lCfgs, E_mMinus1, E_m, nthMin)
            L = 1.1 * L
            continuar = True
        elif factor > 100:
            # choose X with lower energy
            print("C: choosing config with min energy ...")
            # nthMin = 1
            # e, X = select_start_config(lCfgs, E_mMinus1, E_m, nthMin)
            L = 0.9 * L
            continuar = True

        #

        # if ( (I1 < 100) or (I2 < 100)):
        #     assert( len(lCfgs) != 0 )
        #     print("choosing config with min energy ...")
        #     nthMin = 1
        #     e, X = select_start_config(lCfgs, E_mMinus1, E_m, nthMin)
        #     c = counterMin()
        #     if (c % 3 == 2):
        #         print("choosing back ...")
        #         nthMin = 3
        #         e, X = select_start_config(lCfgs, E_mMinus1, E_m, nthMin)
        #     #
        #     # if (c % 5 == 4):
        #     #     print("increasing nSteps ...")
        #     #     nSteps = nSteps * 2
        #     #
        #     continuar = True

        # if (not belongs(ratioAcceptances, ratioMin, ratioMax)):
        #     assert( len(lCfgs) != 0 )
        #     # print("choosing config with max energy ...")
        #     # print("increasing L value ...")
        #     # nthMin = 2
        #     # e, Xtemp = select_start_config(lCfgs, E_mMinus1, E_m, nthMin)
        #     # X = copy.deepcopy(Xtemp)
        #
        #     nthMin = 1
        #     e, X = select_start_config(lCfgs, E_mMinus1, E_m, nthMin)
        #     continuar = True
        #
        #     if ratioAcceptances <= ratioMin:
        #         print("decreasing L value ...")
        #         L = L / 1.2
        #     elif ratioMax <= ratioAcceptances:
        #         print("increasing L value ...")
        #         L = L * 1.2
        #     # continuar = True

    #===========================================================================

    # calculate α_m, log_idos:
    alpha = I_EmToEmplus1 / I_EminToEm # MonteCarlo: ~cociente de goles ∫_Em-1^Em  / ∫_Emin^Em-1
    assert(alpha != 0)
    print( "....................................alpha = " + str(float("{0:.3f}".format(alpha))) + "   L = " + str(float("{0:.3f}".format(L))) )
    log_idos.append( log_sum_idos + log(alpha) ) # =log_idos[m]

    # calculate log_sum_idos:
    # Log of the sum of integrals for each interval, m=1, m=2, m=3,...
    DOSintegrals = [ exp(log_int) for log_int in log_idos ]
    log_sum_idos = log(sum(DOSintegrals))

    return lCfgs, alpha, log_idos, log_sum_idos, L
#
def select_start_config(lCfgs, E_mMinus1, Em, nthMin):
    import copy
    # "During each partitioning process, a set of configurations are randomly
    # stored, from which a correct starting configuration for the next
    # partitioning process is selected.""
    # import copy
    # for i in range(len(lCfgs)):
    #     e, X = lCfgs[i]
    #     if belongs(e, E_mMinus1, Em):
    #         return e, X
    # #
    # #
    # e = lCfgs[0][0]
    # X = copy.deepcopy(lCfgs[0][1])
    # return e, X # I don't know what is the "correct starting configuration"


    # Xsaved = []
    # esaved = []
    # for i in range(len(lCfgs)):
    #     e, X = lCfgs[i]
    #     if belongs(e, E_mMinus1, Em):
    #         esaved.append(e)
    #         Xsaved.append(X)
    # #

    if nthMin == 1:
        Xsaved = []
        esaved = []
        for i in range(len(lCfgs)):
            e, X = lCfgs[i]
            if belongs(e, E_mMinus1, Em):
                esaved.append(e)
                Xsaved.append(X)
        #
        if (len(esaved) > 0):
            minE = min(esaved)
            idx = esaved.index(minE)
            X = Xsaved[idx]
            return minE, X

    elif nthMin == 2:
        # l = sorted(set(esaved))
        # n = int(len(l) / 2.0)
        # n = len(l) - 1
        # minE = l[n] # it is not mininum in fact
        # minE = max(esaved)

        Xsaved = []
        esaved = []
        for i in range(len(lCfgs)):
            e, Xtemp = lCfgs[i]
            X = copy.deepcopy(Xtemp)
            if (Em < e):
                esaved.append(e)
                Xsaved.append(X)
        #
        minE = max(esaved)
        idx = esaved.index(minE)
        Xtemp = Xsaved[idx]
        X = copy.deepcopy(Xtemp)
        return minE, X
    #
    elif nthMin == 3:
        Xsaved = []
        esaved = []
        for i in range(len(lCfgs)):
            e, Xtemp = lCfgs[i]
            X = copy.deepcopy(Xtemp)
            if (e <  E_mMinus1):
                esaved.append(e)
                Xsaved.append(X)
        #
        minE = min(esaved)
        idx = esaved.index(minE)
        Xtemp = Xsaved[idx]
        X = copy.deepcopy(Xtemp)
        return minE, X


    # idx = esaved.index(minE)
    # X = Xsaved[idx]
    # return minE, X


#
#
################################################################################
# comments:
    # log_idos[m] = log_sum_idos + log(alpha)
    # if you develop the sum:
        # = log(  ∫_Emin^Em-1 Ω(E)dE   ) + log(α[m + 1])
        # = log(  ∫_Emin^Em-1 Ω(E)dE *  α[m + 1] )
        # = log(  ∫_Emin^Em-1 Ω(E)dE * ∫_Em-1^Em Ω(E)dE / ∫_Emin^Em-1 Ω(E)dE )
        # = log( ∫_Em-1^Em Ω(E)dE )
        # then, why not just making equal to: log_idos[m] = log(I1) ???
        # ANSWER: because I1 is not equal to ∫_Em-1^Em Ω(E)dE
        # remember, you only know the rate I1/I2 = ∫_Em-1^Em  / ∫_Emin^Em-1
        # the real operation is:
        # = log(  ∫_Emin^Em-1 Ω(E)dE *  α[m + 1] )
        # = log(  ∫_Emin^Em-1 Ω(E)dE * k*∫_Em-1^Em Ω(E)dE / k*∫_Emin^Em-1 Ω(E)dE )
        # where k is an unknown constant!
