import testFunPar as tp
# if rank == 0:
#     for i in range ( size ) :
#         assert data [i] == (i +1) **2
# else :
#     assert data is None
#

data, lista = tp.fun()
print(lista)
# if rank == 0:
myLista = []
for i in range(len(lista)):
    for j in range(len(lista[i])):
        myLista.append(lista[i][j])
    #

print( " [", rank, "] After gather : ", data)
print( " [", rank, "] After gather : ", myLista)
