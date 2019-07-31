from ase import Atom, Atoms
from lammpslib import LAMMPSlib

cmds = ["pair_style eam/alloy",
        "pair_coeff * * NiAlH_jea.eam.alloy Al"]

a = 4.05
al = Atoms([Atom('Al')], cell=(a, a, a), pbc=True)
lammps = LAMMPSlib(lmpcmds = cmds, atom_types={'Al':1}, logfile='test.log')
al.set_calculator(lammps)
print "Energy ", al.get_potential_energy()




pot = LAMMPSlib(lmpcmds=init_cmds, atom_types='H 1', log_file='lammps.%d.log' % rank, keep_alive=True, lammps_name='mpi',
                lammps_header=, lammps_header_extra=["atom_style full"], comm=calculator_comm, read_molecular_info=ns_args['LAMMPS_molecular_info'])



# chinchay@Jaguar in /usr/local/lib/python2.7/site-packages$  ln -s /Users/chinchay/Documents/Instaladores/pymatnest-master/lammpslib.pyc .
from ase import Atom, Atoms
from lammpslib import LAMMPSlib

cmds = ["pair_style eam/alloy",
        "pair_coeff * * /Users/chinchay/Documents/Instaladores/LAMMPS/lammps/potentials/NiAlH_jea.eam.alloy Al H"]

a = 4.05
al = Atoms([Atom('Al')], cell=(a, a, a), pbc=True)
h  = Atoms([Atom('H')])
struct = al + h
struct.get_cell()
struct.get_positions()

# mylammps = LAMMPSlib(lmpcmds = cmds, atom_types={'Al':1, 'H': 1}, logfile='/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/test.log')
mylammps = LAMMPSlib(lmpcmds = cmds, atom_types="TYPE_EQUALS_Z", logfile='/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/test.log')
mylammps.parameters.atom_types
mylammps.parameters.atom_types_equal_atomic_numbers

if isinstance(mylammps.parameters.atom_types,basestring):
   if mylammps.parameters.atom_types == "TYPE_EQUALS_Z":
      mylammps.parameters.atom_types_equal_atomic_numbers = True
      mylammps.parameters.atom_types = {}
      for Z in struct.get_atomic_numbers():
         mylammps.parameters.atom_types[Z] = Z
   else:
      raise ValueError('atom_types parameter "%s" is string, but not TYPE_EQUALS_Z' % mylammps.parameters.atom_types)
else:
   # assume atom_types is a dictionary with symbols (or numbers) as keys
   mylammps.parameters.atom_types_equal_atomic_numbers = False
   symbol_atom_types = mylammps.parameters.atom_types.copy()
   mylammps.parameters.atom_types = {}
   for sym in symbol_atom_types:
      try:
         num = int(sym)
      except:
         num = symbols2numbers(sym)[0]
      mylammps.parameters.atom_types[num] = symbol_atom_types[sym]




[(i,Z) for i,Z in enumerate( struct.get_atomic_numbers() )]
[(i,mylammps.parameters.atom_types[Z]) for i,Z in enumerate( struct.get_atomic_numbers() )]
mylammps.parameters.atom_types






alh.set_calculator(lammps)

alh

print "Energy ", alh.get_potential_energy()


lammps.propagate

alh.get_pbc()
lammps.parameters.atom_types


lammps.propagate(alh, ['energy'], ['positions'], 0)


lammps.results


alh.get_potential_energy()

lmp.extract_variable("eng",None,0)
alh.lmp.extract_variable('pe', None, 0)

print lammps.propagate()
self.lmp.extract_variable('pe', None, 0)

lammps.parameters.create_box

def parseAtomTypes(LAMMPS_atom_types): # LAMMPS_atom_types is a string
    if len(LAMMPS_atom_types) > 0:
        if LAMMPS_atom_types == 'TYPE_EQUALS_Z':
            out = LAMMPS_atom_types
        else:
            out = {}
            for type_pair in [s.strip() for s in LAMMPS_atom_types.split(',')]:
                f = type_pair.split()
                out[f[0]] = int(f[1])
    #
    return out
#
parseAtomTypes('Al 1, H 1')


t = "hola aqui"
t.pop('LAMMPS_atom_types', '')



print(lammps.parameters.atom_types)
lammps.

self.propagate(atoms, properties, system_changes, 0)
self.initialise_lammps(atoms)
self.parameters.atom_types



print "Energy ", alh.get_potential_energy()
# alh.set_calculator(lammps)
# print "Energy ", alh.get_potential_energy()











def main():
        """ Main function """
        global movement_args
        global ns_args, start_first_iter
        global max_n_cull_per_task
        global size, rank, comm, rng, np, sys, ns_analyzers
        global n_cull, n_walkers, n_walkers_per_task
        global n_extra_walk_per_task
        global do_calc_quip, do_calc_lammps, do_calc_internal, do_calc_fortran
        global energy_io, traj_io, walkers
        global n_atoms, KEmax, pot

def propagate_lammps(at, dt, n_steps, algo, Emax=None):
    if pot.first_propagate:
        pot.first_propagate=False
    else:
        pot.lmp.command('unfix 1 ')

    # set appropriate fix
    if algo == 'NVE':
        pot.lmp.command('fix 1 all nve')
    elif algo == 'GMC':
        # hard coded value needed for LAMMPS.  Let's hope our RNG maximum is at least this large
        pot.lmp.command('fix 1 all gmc {} {} '.format(rng.int_uniform(1,900000000),Emax))
    else:
        exit_error("propagate_lammps got algo '%s', neither NVE nor GMC\n"% algo)

    # NB maybe not do this every time? Just _after_ MD, since that's the only way position can become crazy?
    at.wrap()

    # do propagation, recover from exception if needed
    try:
        if algo == 'NVE':
            pot.propagate(at, properties=['energy','forces'],system_changes=['positions'], n_steps=n_steps, dt=dt)
        else:
            pot.propagate(at, properties=['energy','forces'],system_changes=['positions'], n_steps=n_steps, dt=dt, dt_not_real_time=True)
    except Exception as err:
        # clean up and return failure
        if ns_args['debug'] >= 4:
            print "propagate_lammps got exception ", err
        pot.restart_lammps(at)
        pot.first_propagate=True
        return False

    return True
#

            if propagate_lammps(at, dt=movement_args['MD_atom_timestep'], n_steps=movement_args['atom_traj_len'], algo='NVE'):
                final_E = pot.results['energy'] + eval_energy(at, do_PE=False)
            else: # propagate returned success == False
                final_E = 2.0*abs(Emax)
                ## print "error in propagate_lammps NVE, setting final_E = 2*abs(Emax) =" , final_E





        elif ns_args['energy_calculator'] == 'lammps':
            try:
                from lammpslib import LAMMPSlib
            except:
                exit_error("energy_calculator=lammps and failed to import lammpslib module\n", 1)
            do_calc_lammps=True
            try:
                ns_args['LAMMPS_init_cmds'] = args.pop('LAMMPS_init_cmds')
            except:
                exit_error("need LAMMPS initialization commands LAMMPS_init_cmds\n",1)
            ns_args['LAMMPS_name'] = args.pop('LAMMPS_name', '')
            ns_args['LAMMPS_header'] = args.pop('LAMMPS_header', 'units metal; atom_style atomic; atom_modify map array sort 0 0')
            ns_args['LAMMPS_header_extra'] = args.pop('LAMMPS_header_extra', '')
            ns_args['LAMMPS_atom_types'] = None

            ns_args['LAMMPS_fix_gmc'] = str_to_logical(args.pop('LAMMPS_fix_gmc', "F"))
            LAMMPS_atom_types = args.pop('LAMMPS_atom_types', '')
            if len(LAMMPS_atom_types) > 0:
                if LAMMPS_atom_types == 'TYPE_EQUALS_Z':
                    ns_args['LAMMPS_atom_types'] = LAMMPS_atom_types
                else:
                    ns_args['LAMMPS_atom_types'] = {}
                    for type_pair in [s.strip() for s in LAMMPS_atom_types.split(',')]:
                        f = type_pair.split()
                        ns_args['LAMMPS_atom_types'][f[0]] = int(f[1])
            else:
               exit_error("LAMMPS_atom_types is mandatory if calculator type is LAMMPS\n",1)





        # initialize mpi
        comm = None
        calculator_comm = None
        rank = 0
        size = 1
        if use_mpi:
            print "INFO: use_mpi true, importing mpi4py module"
            try:
                from mpi4py import MPI
            except:
                sys.stderr.write("Failed to import mpi4py\n")
                sys.exit(10)
            comm = MPI.COMM_WORLD
            calculator_comm = MPI.COMM_SELF

        if use_mpi:
            try:
                rank = comm.Get_rank()
                size = comm.Get_size()
            except:
                exit_error("Failed to get rank or size\n", 10)

        if comm is not None:
            print "comm ", comm, " size ", size, " rank ", rank

        # read inputs on root, then bcast
        if rank == 0:
            lines=sys.stdin.readlines()
            if len(lines) == 0:
                try:
                    infile=open("ns_inputs","r")
                except:
                    exit_error("Failed to read ns_inputs file\n", 1)
                lines = infile.readlines()
            args={}
            if rank == 0:
                for line in lines:
                    if re.match("\s*(#.*)?$", line):
                        continue
                    matches = re.match("\s*(\S+)\s*=\s*(.*\S)", line)
                    if matches is None:
                        exit_error("Failed to parse line '%s'" % line, 1)
                    args[matches.group(1)] = matches.group(2)
        else:
            args = None
        if comm is not None:
            args = comm.bcast(args,root=0)


init_cmds = [s.strip() for s in ns_args['LAMMPS_init_cmds'].split(';')]
header_cmds = [s.strip() for s in ns_args['LAMMPS_header'].split(';')]

pot = LAMMPSlib(lmpcmds=init_cmds, atom_types=ns_args['LAMMPS_atom_types'], log_file='lammps.%d.log' % rank, keep_alive=True, lammps_name=ns_args['LAMMPS_name'],
                lammps_header=header_cmds, lammps_header_extra=header_extra_cmds, comm=calculator_comm, read_molecular_info=ns_args['LAMMPS_molecular_info'])


pot = LAMMPSlib(lmpcmds=init_cmds, atom_types='H 1', log_file='lammps.%d.log' % rank, keep_alive=True, lammps_name='mpi',
                lammps_header=header_cmds, lammps_header_extra=header_extra_cmds, comm=calculator_comm, read_molecular_info=ns_args['LAMMPS_molecular_info'])


print "PRE START_LAMMPS"
pot.start_lammps() # so that top level things like units will be set
print "POST START_LAMMPS"
pot.first_propagate=True


at.set_calculator(pot)




if not ns_args['LAMMPS_atom_types'] == 'TYPE_EQUALS_Z':
    try:
        used_chem_symbols = {ase.atoms.generalized_chemical_symbols(int(species.split()[0])) for species in species_list}
    except:
        used_chem_symbols = {ase.atoms.chemical_symbols[int(species.split()[0])] for species in species_list}
    if not used_chem_symbols == set(ns_args['LAMMPS_atom_types'].keys()):
        exit_error("species in start_species must correspond to those in LAMMPS_atom_types\n",1)






if rank == 0:
    energy_io.close()
traj_io.close()
if track_traj_io is not None:
    track_traj_io.close()
if E_dump_io is not None:
    E_dump_io.close()

if comm is not None:
    MPI.Finalize()
sys.exit(0)
