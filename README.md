# lavax++
LAmmps VAsp eXchanger, written in C++.

Dependencies for compilation are <code>boost</code> and <code>blitz++</code>.

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

## lavax.conf
lavax must be setup according to the following configuration file options:

```squidconf
VASP_COMMAND    = mpirun -np 32 vasp533_mpi
LAMMPS_COMMAND  = mpirun -np 8 lmp_mpi
INIT_POSCAR     = POSCAR_init         # File describing the initial crystal state
LVX_ITERATIONS  = 70                  # Number of internal LAVAX iterations.
LAMMPS_POTENTIAL_FILE = W_BN.eam.fs

# Options for adaptive timestep in VASP and LAMMPS:
USE_ADAPTIVE_TIMESTEP = true
MAX_DISTANCE = 0.1 # Angstroms
MAX_TIMESTEP = 3   # Femtoseconds

# The distance where the good and bad potentials starts to depart:
POTENTIAL_DEPARTURE_DISTANCE = 2.0

# List all the atomic elements in the simulation cell:
ATOMIC_SYMBOL_0 = W
ATOMIC_MASS_0   = 184.0
VASP_POTENTIAL_FILE_GOOD_0   = POTCAR0_good # POTCAR will be concatenated
VASP_POTENTIAL_FILE_BAD_0    = POTCAR0_bad  # from these.
VASP_POTENTIAL_SYMBOL_GOOD_0 = Wsv
VASP_POTENTIAL_SYMBOL_BAD_0  = W
```

## POSCAR_init
The initial POSCAR must contain a line containing the potential symbols, specified in the lavax.conf file.
```text
BCC Xx 
1.0
6.36000000 0.00000000 0.00000000
0.00000000 6.36000000 0.00000000
0.00000000 0.00000000 6.36000000
W Wsv                             <- This line must not be omitted!
14 2
Cartesian
0.00000000 0.00000000 0.00000000
1.59000000 1.59000000 1.59000000
...
```
