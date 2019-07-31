def fun():
    from mpi4py import MPI
    comm  = MPI . COMM_WORLD
    size  = comm . size
    rank  = comm . rank
    data  = ( rank +1) **2
    lista = [data, data + 1]
    print( " [", rank, "] Sending value : ", data)

    data  = comm.gather( data , root =0)
    lista = comm.gather( lista, root =0)
    return data, lista
