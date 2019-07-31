
def initializeParameters(aseStruct, path, pot):
    # path should not end in '/'
    # make these variables global:
    epsilon = 0.1
    sigma   = 2.5
    rcut    = 7.50
    aCell   = 15.0
    # path = '/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning'
    lammps_header = [ "units       metal"]
    cmdsLJ        = [ "pair_style  mlip  " + path + "/mlip_LJ.ini" ,\
                      "pair_coeff  * * " ]
    cmdsMTP       = [ "pair_style  mlip  " + path + "/mlip.ini" ,\
                      "pair_coeff  * * " ]
    log_file      = path + "/log.txt"
    # defining dictionary of Z:LAMMPStype
    atom_types = {}
    for i in range(len(listZ)):
        atom_types[listZ[i]] = i + 1
    #
    if pot == 'mtp':
        cmds = [ "pair_style  mlip  " + path + "/mlip.ini" ,\
                 "pair_coeff  * * " ]
    elif pot == 'aseLJ':
        cmds = [ "pair_style  mlip  " + path + "/mlip_LJ.ini" ,\
                 "pair_coeff  * * " ]
    #
    return atom_types, epsilon, sigma, rcut, aCell, path, lammps_header, cmds, log_file
#
# make these variables global:
atom_types, epsilon, sigma, rcut, aCell, path, lammps_header, cmds, log_file = initializeParameters(aseStruct, path, pot)


def getEnergyAndEigen(aseStruct, isForceMatrixRequired):
    # the following parameteres were defined when initializeParameters() was called.
    mylammps = LAMMPSlib(lmpcmds=cmds, atom_types=atom_types, keep_alive=True, log_file=log_file)
    aseStruct.set_calculator(mylammps)
    energy = aseStruct.get_potential_energy()
    eig = []
    if isForceMatrixRequired:
        ph = Phonons(aseStruct, mylammps)
        ph.run()
        ph.read(acoustic=True)
        ph.clean()
        f = ph.get_force_constant()
        (l,m,n) = f.shape
        if l == 1:
            ff = np.reshape(f, (m,n))
        else:
             print("error")
        #
        eig = LA.eigvalsh( ff ) # eig is a numpy array
    #
    return energy, [float("{0:.5f}".format(eig[i])) for i in range(len(eig))]
#
