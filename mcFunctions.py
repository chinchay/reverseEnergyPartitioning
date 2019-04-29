#call mc_move():
    #get Xold -> Xnew, small enough such that the `E` resultant does not go outside of [Emin,Em]
    #accept the movement Xold->Xnew with probability P(o->n)=min(1, w(En)/w(Eo))
    #get the new energy of Xnew: E

from math import log, exp
# log2 = log(2)

from functions import belongs
from matplotlib.pylab import *  # for random()

def getEm(α_mMinus1, E_mMinus1, E_mMinus2):
    # minimo = min( [ α[m - 1], 1 ] )
    # rLog   = log( 1 +  (1 / minimo) ) / log2
    # Em     = E[m - 1] + (  (E[m - 1] - E[m - 2]) * rLog  )
    minimo = min([ 1, α_mMinus1 ])
    # log3 = log2 + 1.0
    rLog   = log( 1 +  (1 / minimo) ) / log(2) #<<<<<<<<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<<<<<<<
    return E_mMinus1 + (  (E_mMinus1 - E_mMinus2) * rLog  )
#

def getRandomWalk(L): # L: amplitud
    return L * (random() - 0.5) # `0.5`:for centering around 0.

def getRandomWalk3D(L):
    # L = 0.1 # think a better factor <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!!
    # L is related to ~diffusion coefficient??? L must be small enough
    return [getRandomWalk(L), getRandomWalk(L), getRandomWalk(L)] # = [dx,dy,dz]
    # return [getRandomWalk(L), getRandomWalk(L), 0] # = [dx,dy,dz] <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#

def addVectors(a, b): # `a` and `b` must have the same lenght!
    from vectormath import Vector3 as Vector3
    return Vector3([a[i] + b[i] for i in range(len(a))])


def getNewConfig(i, X, L): # L ~amplitud of random walk
    # X_old is a vector with [r1,r2,r3,...,rNatoms] with the positions of every
    # atom of the configuration. ri = [x,y,z]
    iAtom = i % len(X)

    # # Apply periodicity for deltaPosition:
    # if X_old[iAtom].isSlavePeriodic:
    #     jMaster = X_old[iAtom].indexMaster
    #     X_old[iAtom].position = X_old[jMaster].position
    #     return X_old
    # meanwhile, instead of passing atoms:

    # if (iAtom == 5):
    #     dX = getRandomWalk3D(L)
    #     X_old[iAtom] = addVectors( X_old[iAtom], dX )
    #     X_old[9]    = addVectors( X_old[9],    dX )
    #     return X_old
    # if (iAtom == 9):
    #     return X_old
    #
    # if (iAtom == 10):
    #     dX = getRandomWalk3D(L)
    #     X_old[iAtom] = addVectors( X_old[iAtom], dX )
    #     X_old[14]    = addVectors( X_old[14],    dX )
    #     return X_old
    # if (iAtom == 14):
    #     return X_old
    #
    # if (iAtom == 15):
    #     dX = getRandomWalk3D(L)
    #     X_old[iAtom] = addVectors( X_old[iAtom], dX )
    #     X_old[19]    = addVectors( X_old[19],    dX )
    #     return X_old
    # if (iAtom == 19):
    #     return X_old
    # ####
    #
    #
    # if (iAtom == 1):
    #     dX = getRandomWalk3D(L)
    #     X_old[iAtom] = addVectors( X_old[iAtom], dX )
    #     X_old[21]    = addVectors( X_old[21],    dX )
    #     return X_old
    # if (iAtom == 21):
    #     return X_old
    #
    # if (iAtom == 2):
    #     dX = getRandomWalk3D(L)
    #     X_old[iAtom] = addVectors( X_old[iAtom], dX )
    #     X_old[22]    = addVectors( X_old[22],    dX )
    #     return X_old
    # if (iAtom == 22):
    #     return X_old
    #
    # if (iAtom == 3):
    #     dX = getRandomWalk3D(L)
    #     X_old[iAtom] = addVectors( X_old[iAtom], dX )
    #     X_old[23]    = addVectors( X_old[23],    dX )
    #     return X_old
    # if (iAtom == 23):
    #     return X_old
    #
    # #####
    #
    # if (iAtom == 0):
    #     dX = getRandomWalk3D(L)
    #     X_old[iAtom] = addVectors( X_old[iAtom], dX )
    #     X_old[4]     = addVectors( X_old[4],    dX )
    #     X_old[20]    = addVectors( X_old[20],    dX )
    #     X_old[24]    = addVectors( X_old[24],    dX )
    #     return X_old
    # if ((iAtom == 4) or (iAtom == 20) or (iAtom == 24)):
    #     return X_old


    X[iAtom] = addVectors( X[iAtom], getRandomWalk3D(L) )
    return X
#

################################################################################
def getEnergyConfig(X, a1, a2, a3):
    import potential
    from potential import VLJ as V
    import harmonic as ha
    # E = 0

    E, energyPerParticle = ha.get3DPotNearFirstNeighb(X, V, a1, a2, a3)

        # E = potential.getTotalPotentialEnergy(X, V2)
    # E = potential.get3DPotNearFirstNeighb(X, V, a1, a2, a3)
        # E = 0.2*(random()-0.5) - 0.56  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # print(E)
        # E = getEnergyFromStaticMTP(X) # do not relax!  <<<<<<<<<<<<<<<<<<<<<<<<<<<
    return E, energyPerParticle
#
################################################################################


def getWeigth(E, E_m_minus_1, E_m_minus_2):
    # see the correct formula in https://pubs.acs.org/doi/ipdf/10.1021/ct3007056
    if (E_m_minus_1 <= E):
        f = -3 * log(2)
        g = (E - E_m_minus_1) / (E_m_minus_1 - E_m_minus_2)
        return exp(f * g)
    else: # (E < E_m_minus_1)
        return 1.0
#

def getProbTransition(E_new, E_old, E_m_minus_1, E_m_minus_2):
    w_new = getWeigth(E_new, E_m_minus_1, E_m_minus_2)
    w_old = getWeigth(E_old, E_m_minus_1, E_m_minus_2)
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

def mc_move(i, E_old, X_old, Xeq, forceMatrix, E_m_minus_1, E_m_minus_2, L, a1, a2, a3):
    import copy
    import harmonic as ha



    Xtemp  = copy.deepcopy(X_old)
    X_temp = getNewConfig(i, Xtemp, L)
    X_new  = copy.deepcopy(X_temp)
    E_new, energyPerParticle = getEnergyConfig(X_new, a1, a2, a3)
    P_old2new = getProbTransition(E_new, E_old, E_m_minus_1, E_m_minus_2)
    hasMoved = False

    if ( P_old2new >= random() ):
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
            rHyper, deltaEharmonic = ha.getHarmonicEnergy(X_new, Xeq, forceMatrix)
            out = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/rHyperEHarmonic" + str(c) + ".txt"
            s = str(rHyper) + " , " + str(deltaEharmonic) + " , " +  str(E_new) + "\n"
            # s = str(rHyper) + " , " + str(deltaEharmonic) + "\n"
            f = open(out, "w")
            f.write(s)
            f.close()


        return E_new, X_new, hasMoved
    #
    X = copy.deepcopy(X_old)
    return E_old, X, hasMoved

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
def randomMCmoves(Eminimo, E_mMinus2, E_mMinus1, E_m,\
                  nSteps, e, X, log_idos, log_sum_idos, L, a1, a2, a3, Xeq, forceMatrix):
    import copy
    import harmonic as ha


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
    repMax = 1 #10 # maximum number of repetitions to get the appropriate ratio of
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

        # if (repeat > 0):
        #     if I2 == 0:
        #         repeat = 0
        #         # L = 0.99 * L
        #     # elif ratioAcceptances > ratioMax: L = 1.01 * L
        #     # elif ratioAcceptances < ratioMax: L = 0.99 * L


        [I1, I2] = [0, 0]
        repeat += 1
        acceptances = 0

        # L = L * 0.8 # Decrease the amplitud of the random walk to get the
                    # desired ratioAcceptances range values [25%, 35%].
                    # You can improve with a more sophisticated algorithm here.

        # Collect ehist for [Emin,E[m-1]] and [E[m-1],Em]:
        for i in range(nSteps):
            # get a configuration e,X after randomly move X:

            e, X, hasMoved = mc_move(i, e, X, Xeq, forceMatrix, E_mMinus1, E_mMinus2, L, a1, a2, a3)

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
                if belongs(e, E_mMinus1, E_m): # (E[m - 1] <= e <= Em):
                    I1 += 1  # un gol más para ∫_Em-1^Em Ω(E)dE
                elif belongs(e, Eminimo, E_mMinus1): # (Emin <= e <= E[m - 1]): #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                # elif belongs(e, E_mMinus1 - 10, E_mMinus1): # (Emin <= e <= E[m - 1]):
                    I2 += 1  # un gol más para ∫_Emin^Em-1 Ω(E)dE
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
        print(ratioAcceptances, I2, I1, e, E_mMinus1, E_m, belongs(e, Eminimo, E_mMinus1), belongs(e,E_mMinus1, E_m), L, repeat)


        factor =  (I1 + 1) / (I2 + 1) # '+1' to avoid zeros
        assert( len(lCfgs) != 0 )
        if ( (I1 < 100) or (I2 < 100)):
            # choose X with lower energy
            print("A: choosing config with min energy ...")
            nthMin = 1
            e, X = select_start_config(lCfgs, E_mMinus1, E_m, nthMin)
            continuar = True
        elif factor < 1/100.0:
            # choose X with higher energy
            print("B: choosing config with high energy ...")
            nthMin = 2
            e, X = select_start_config(lCfgs, E_mMinus1, E_m, nthMin)
            continuar = True
        elif factor > 100:
            # choose X with lower energy
            print("C: choosing config with min energy ...")
            nthMin = 1
            e, X = select_start_config(lCfgs, E_mMinus1, E_m, nthMin)
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
    alpha = I1 / I2 # MonteCarlo: ~cociente de goles ∫_Em-1^Em  / ∫_Emin^Em-1
    assert(alpha != 0)
    print(alpha, L)
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
