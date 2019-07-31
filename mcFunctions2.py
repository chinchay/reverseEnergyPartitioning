################################################################################
# python3
import os
print(os.getcwd())
os.chdir('/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/')
print(os.getcwd())

# import graphlab
################################################################################

#call mc_move():
    #get Xold -> Xnew, small enough such that the `E` resultant does not go outside of [Emin,Em]
    #accept the movement Xold->Xnew with probability P(o->n)=min(1, w(En)/w(Eo))
    #get the new energy of Xnew: E

from math import log, exp
# log2 = log(2)

import numpy as np

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

# def addVectors(a, b): # `a` and `b` must have the same lenght!
#     from vectormath import Vector3 as Vector3
#     return Vector3([a[i] + b[i] for i in range(len(a))])


def getNewConfig(i, structX, L): # L ~amplitud of random walk
    # X_old is a vector with [r1,r2,r3,...,rNatoms] with the positions of every
    # atom of the configuration. ri = [x,y,z]
    # NOW IT IS IN FRACTIONAL COORDINATES!! https://pymatgen.org/ "# Changes species and coordinates (fractional assumed for structures)" CARTESIANS for molecules
    iAtom = i % structX.num_sites
    structX[iAtom] = structX[iAtom].specie, structX[iAtom].frac_coords + np.array(getRandomWalk3D(L))
    return structX
#

################################################################################
def getEnergyConfig(structX, radio, dRneighborsFromEachSite):
    import potential
    from potential import VLJ as V
    # import harmonic as ha
    import crystal as cr
    # E = 0

    # E, energyPerParticle = ha.get3DPotNearFirstNeighb(structX, V, a1, a2, a3)
    energy = cr.getEnergy(structX, radio, dRneighborsFromEachSite)

        # E = potential.getTotalPotentialEnergy(X, V2)
    # E = potential.get3DPotNearFirstNeighb(X, V, a1, a2, a3)
        # E = 0.2*(random()-0.5) - 0.56  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # print(E)
        # E = getEnergyFromStaticMTP(X) # do not relax!  <<<<<<<<<<<<<<<<<<<<<<<<<<<
    return energy
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

def mc_move(i, E_old, structX_old, structXeq, forceMatrix, E_m_minus_1, E_m_minus_2, L, radio, dRneighborsFromEachSite):
    import copy
    import harmonic as ha

    structXtemp  = copy.deepcopy(structX_old)
    structX_temp = getNewConfig(i, structXtemp, L)
    structX_new  = copy.deepcopy(structX_temp)
    E_new  = getEnergyConfig(structX_new, radio, dRneighborsFromEachSite)
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
            nAtoms = structX_old.num_sites
            out = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/moves." + str(c) + ".xyz"
            pos = str(nAtoms) + "\n\n"
            for i in range(nAtoms):
                pos += str(structX_new[i].specie) + str(structX_new[i].x) + "   " + str(structX_new[i].y) + "   " + str(structX_new[i].z) + "\n"
            #
            pos += "\n"

            f = open(out, "w")
            f.write(pos)
            f.close()

            #######################################################
            rHyper, deltaEharmonic = ha.getHarmonicEnergy(structX_new.cart_coords, structXeq.cart_coords, forceMatrix)
            out = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/rHyperEHarmonic" + str(c) + ".txt"
            s = str(rHyper) + " , " + str(deltaEharmonic) + " , " +  str(E_new) + "\n"
            # s = str(rHyper) + " , " + str(deltaEharmonic) + "\n"
            f = open(out, "w")
            f.write(s)
            f.close()


        return E_new, structX_new, hasMoved
    #
    structX = copy.deepcopy(structX_old)
    return E_old, structX, hasMoved

################################################################################
# MC sampling with weight w(E). Collect quantities for the next subdivision m+1.
################################################################################
def randomMCmoves(Eminimo, E_mMinus2, E_mMinus1, E_m,\
                  nSteps, e, structX, log_idos, log_sum_idos, L, structXeq, forceMatrix, radio, dRneighborsFromEachSite):
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

            e, structX, hasMoved = mc_move(i, e, structX, structXeq, forceMatrix, E_mMinus1, E_mMinus2, L, radio, dRneighborsFromEachSite)

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
                if random() > ratioOfCfgsToSave: lCfgs.append([e, structX])
            #
            if (e < 1.1 * Eminimo):
                print(e, Eminimo)
                ########################################################################################
                import os
                c = counter()
                nAtoms = structX_old.num_sites
                out = "/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/ovitoPlots/moves." + str(c) + ".xyz"
                pos = str(nAtoms) + "\n\n"
                for i in range(nAtoms):
                    pos += str(structX_new[i].specie) + str(structX_new[i].x) + "   " + str(structX_new[i].y) + "   " + str(structX_new[i].z) + "\n"
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
            e, structX = select_start_config(lCfgs, E_mMinus1, E_m, nthMin)
            continuar = True
        elif factor < 1/100.0:
            # choose X with higher energy
            print("B: choosing config with high energy ...")
            nthMin = 2
            e, structX = select_start_config(lCfgs, E_mMinus1, E_m, nthMin)
            continuar = True
        elif factor > 100:
            # choose X with lower energy
            print("C: choosing config with min energy ...")
            nthMin = 1
            e, structX = select_start_config(lCfgs, E_mMinus1, E_m, nthMin)
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

    if nthMin == 1:
        structXsaved = []
        esaved = []
        for i in range(len(lCfgs)):
            e, structX = lCfgs[i]
            if belongs(e, E_mMinus1, Em):
                esaved.append(e)
                structXsaved.append(structX)
        #
        if (len(esaved) > 0):
            minE = min(esaved)
            idx = esaved.index(minE)
            structX = structXsaved[idx]
            return minE, structX

    elif nthMin == 2:
        # l = sorted(set(esaved))
        # n = int(len(l) / 2.0)
        # n = len(l) - 1
        # minE = l[n] # it is not mininum in fact
        # minE = max(esaved)

        structXsaved = []
        esaved = []
        for i in range(len(lCfgs)):
            e, structXtemp = lCfgs[i]
            X = copy.deepcopy(structXtemp)
            if (Em < e):
                esaved.append(e)
                structXsaved.append(structX)
        #
        minE = max(esaved)
        idx = esaved.index(minE)
        structXtemp = structXsaved[idx]
        structX = copy.deepcopy(structXtemp)
        return minE, structX
    #
    elif nthMin == 3:
        structXsaved = []
        esaved = []
        for i in range(len(lCfgs)):
            e, structXtemp = lCfgs[i]
            structX = copy.deepcopy(structXtemp)
            if (e <  E_mMinus1):
                esaved.append(e)
                structXsaved.append(structX)
        #
        minE = min(esaved)
        idx = esaved.index(minE)
        structXtemp = structXsaved[idx]
        structX = copy.deepcopy(structXtemp)
        return minE, structX


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
