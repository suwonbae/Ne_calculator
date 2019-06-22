# NePP

Calculates Ne from a primitive path network even with cross-links considered.

## Getting Started

### dump.PP
A dump file corresponding a primitive path network, which is an input file to be read in.

### run.sh
Runs the main program (Neq_mpi.f) on every single primitive path.
User changes two variables, num_chains and num_procs corresponding to the total number of macromolecular chains and processors.

### Neq_mpi.f
The main program that generates a fine mesh (3D grid) for every single primitive path and counts the number of grid points on which primitive path falls in order to not double count some segments.

### var.default
Default parameters.

### var.options
Parameters that user changes.

#### num_atoms
The total number of atoms in system.

#### num_atoms_per_mol
The total number of atoms in one molecule.

#### griddel
The increment to generates grid points.

#### rprobe
The probe radius describing the radius of a tube confining a macromolecular chain. Only grid points falling within this radius from primitive path are counted.

### postprocess.py
Ne is a statistical property, so Ne is calculated by using the contour length of all primitive paths.
User changes two variables, num_chains and num_atoms_per_mol corresponding to the total number of macromolecular chains and atoms in one molecule.

## Running
```
mpif90 Neq_mpi.f
```
This compiles the program wrttien in fortran and generates a.out.

```
./run.sh
```
This generates stdout's as many as the number of macromolcular chains.

```
python postprocess.py
```
This calculates a statistical property, Ne, by averaging over stdout's.

##
@article{doi:10.1021/acs.macromol.8b01027,	
author = {Bae, Suwon and Galant, Or and Diesendruck, Charles E. and Silberstein, Meredith N.},	
title = {The Effect of Intrachain Cross-Linking on the Thermomechanical Behavior of Bulk Polymers Assembled Solely from Single Chain Polymer Nanoparticles},	
journal = {Macromolecules},	
volume = {51},	
number = {18},	
pages = {7160-7168},	
year = {2018},
doi = {10.1021/acs.macromol.8b01027},
URL = {
https://doi.org/10.1021/acs.macromol.8b01027
}

}
