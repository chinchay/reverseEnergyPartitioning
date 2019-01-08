# Comment

import potential

lista = potential.getAtoms()
a = lista[0]
b = lista[1]
c = lista[2]
a.position
b.position
c.position
v = b.position - c.position
v.length

potential.VLJ(b,c)

potential.test()

def fun(x::Int):
    return 0
#
