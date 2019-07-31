from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print("hello from ", str(rank), " in ", str(size) )

def

#### Parallel work: ########################################################
I1, I2, lCfgs = mc.randomMCmovesParallel2(Emin, E[m - 2], E[m - 1], E[m], nSteps, e, Xcopy, L, structEq, Xeq, forceMatrix)
I1 = comm.gather(I1, root =0)
I2 = comm.gather(I2, root =0)
lCfgs = comm.gather(lCfgs, root =0)
print("end of first Parallel")
############################################################################
