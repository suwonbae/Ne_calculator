# NePP

Calculates Ne from a primitive path network even with cross-links considered.

## Getting Started

### dump.PP
A dump file corresponding a primitive path network, which is an input file to be read in.

### run.sh
Runs the main program (Neq_mpi.f) on every single primitive path.
User changes two variables, num_chains and num_procs corresponding to the total number of macromolecular chain and processors.

### Neq_mpi.f
The main program that generates a fine mesh (3D grid) for every single primitive path and counts the number of grid points on which primitive path falls in order to not double count some segments.

### inc.default
Default parameters.

### inc.options
Parameters that user changes.

#### num_atoms
The total number of atoms in system.

#### num_atoms_per_mol
The total number of atoms in one molecule.

#### griddel
The increment to generates grid points.

#### rvdw
The probe radius describing the radius of tube confining macromolecular chain. Only grid points falling within this radius from primitive path are counted.

### postprocess.py
Ne is a statistical property, so Ne is calculated by using the contour length of all primitive paths.
User changes two variables, num_chains and num_atoms_per_mol corresponding to the total number of macromolecular chains and atoms in one molecule.

## Running
```
./run.sh
```
This generates stdout's as many as the number of macromolcular chains.

```
python postprocess.py
```
This calculates a statistical property, Ne, by averaging over stdout's.
