# POTENTIAL

import Pkg
Pkg.add("LinearAlgebra")
using LinearAlgebra # for norm()


function VLJ(r)
    """
    Lennard-Jones potential
    ε  is the depth of the potential well.
    At rm, the potential function has the value −ε.
    """
    rm = 1.0
    ε  = 1.0

    # Pauli repulsion at short ranges due to overlapping electron orbitals
    repulsive  = (rm / r) ^ 12

    # attraction at long ranges (van der Waals force, or dispersion force??)
    attractive = -2 * ( (rm / r) ^ 6 )
    return ε * ( repulsive + attractive )
end

# abstract type Point end


module Atom
    function getX()
        return 0
    end
end

mutable struct Atomo
    name :: String
    x :: Float16
    y :: Float16
    z :: Float16
    # still I cannot find how to define setter, getters, methods,
    # non-default constructors, inheritance!!
end

function getAtomsDist(A::Atomo, B::Atomo)
    return norm( [A.x - B.x, A.y - B.y, A.z - B.z ] )
end

function atoms()
    A = Atomo("Co", 0.0, 0.0, 0.0)
    B = Atomo("Co", 0.0, 1.0, 0.0)
    vAtoms = [A,B]

    VLJ(r)


end



#
# function modelParticles(N)
#     for i in 1:N
