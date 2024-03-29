# lavax++
LAmmps VAsp eXchanger, written in C++.

Dependencies for compilation are <code>boost</code>.

To compile lavax first run
```bash
autoreconf -i
```
Configure the environment
```bash
./configure --prefix=/desired/install/directory
```
then
```bash
make
```
```bash
make install
```
On some systems it might be necessary to include compiler flags for non-standard installation locations for boost:
```bash
make CXXFLAGS='-I/cfs/klemming/scratch/d/danielk5/local/include -L/cfs/klemming/scratch/d/danielk5/local/lib'
```
Alternatively, you can set the following shell variables instead:
```bash
export CPATH=$CPATH:/cfs/klemming/scratch/d/danielk5/local/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cfs/klemming/scratch/d/danielk5/local/lib
export LIBRARY_PATH=$LIBRARY_PATH:/cfs/klemming/scratch/d/danielk5/local/lib
```
Make the manual via
```bash
make lavax-manual
```
for a comprehensive guide.

## TODO
  * ~~Implement LAMMPS script for multiple atomic species.~~
  * ~~Test adaptive number of NSW steps.~~
  * Remove energy from system at simulation cell boundary.

## lavax.conf
lavax must be setup according to the following configuration file options:

```squidconf
VASP_COMMAND    = mpirun -np 32 vasp533_mpi
LAMMPS_COMMAND  = mpirun -np 8 lmp_mpi
INIT_POSCAR     = POSCAR_init         # File describing the initial crystal state.
LAVAX_ITERATIONS  = 70                # Number of internal LAVAX iterations.
LAMMPS_POTENTIAL_FILE = W_BN.eam.fs

LAMMPS_HIDE_OUTPUT = false
VASP_HIDE_OUTPUT   = false
REFORMAT_XDATCAR   = true

# Options for adaptive timestep in VASP and LAMMPS:
USE_ADAPTIVE_POTIM = true
MAX_DISTANCE = 0.1 # Angstroms
MAX_TIMESTEP = 3   # Femtoseconds

# Switching too many potentials may cause VASP to diverge:
MAX_POTENTIAL_SUBSTITUTIONS = 4
MAX_VASP_NSW = 10  # Maximum number of VASP timesteps before restart

# The distance where the hard and soft potentials starts to depart:
POTENTIAL_DEPARTURE_DISTANCE = 2.0

# List all the atomic elements in the simulation cell:
ATOMIC_SYMBOL_0 = W
VASP_POTENTIAL_FILE_HARD_0   = POTCAR0_hard # POTCAR will be concatenated
VASP_POTENTIAL_FILE_SOFT_0   = POTCAR0_soft # from these.
VASP_POTENTIAL_SYMBOL_HARD_0 = Ws
VASP_POTENTIAL_SYMBOL_SOFT_0 = W
LAMMPS_ATOMIC_SYMBOL_0       = W
```

## POSCAR_init
The initial POSCAR must contain a line with the potential symbols, as specified in the lavax.conf file.
```text
BCC Xx 
1.0
6.36000000 0.00000000 0.00000000
0.00000000 6.36000000 0.00000000
0.00000000 0.00000000 6.36000000
W Ws                              <- This line must not be omitted!
14 2
Cartesian
0.00000000 0.00000000 0.00000000
1.59000000 1.59000000 1.59000000
...
```
## Caveats
The symbols for the atomic elements (VASP_POTENTIAL_SYMBOL_*) should not be more than two characters long due to a bug in VASP.
