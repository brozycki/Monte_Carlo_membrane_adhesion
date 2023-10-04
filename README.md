# Monte Carlo code to simulate membrane adhesion

This C++ code can be used to perform Monte Carlo simulations of the lattice-based mesoscale model introduced in:
- Long Li, Jinglei Hu, Bartosz Rozycki, Fan Song, Nano Letters 20(1):722-728 (2020)
- Long Li, Jinglei Hu, Xinghua Shi, Bartosz Rozycki, Fan Song, Soft Matter 17(7):1912-1920 (2021)

Authors: Lukasz Milewski & Bartosz Rozycki, Institute of Physics, Polish Academy of Sciences, October 2023

This code supplements a book chapter titled "Lattice-based mesoscale simulations and mean-field theory of cell membrane adhesion" by Long Li, Lukasz Milewski, Jie Gao, Jinglei Hu and Bartosz Rozycki. Book title: "Biophysical approaches to lateral and transverse lipid membrane heterogeneity". Editors: Markus Deserno and Tobias Baumgart. 

The development of this Monte Carlo simulation code was supported by the National Science Centre of Poland via grant number 2021/40/Q/NZ1/00017.

# Start

Compile the C++ code 

c++ MC_membrane_adhesion.cpp -O -o sim.o

and execute it with 11 parameters at input as in the following example:

./sim.o 0.02 0.1 10.0 1.5 3.0 6.0 10000000 100000000 1000 20000000 1

The input parameters needs to be given in the following order:

1: area concentration of protein particles (in units of 1/a^2)

2: fraction of membrane area occupied by lipid rafts

3: membrane bending rigidity modulus (in kT units)

4: energy of contact between rafts (in kT units)

5: energy of association of a protein particle with a raft (in kT units)

6: energy of receptor-ligand binding (in kT units)

7: number of MC cycles for equilibration

8: number of MC cycles used for data acquisition

9: number of MC cycles between recording observables

10: number of MC cycles between recording system configurations

11: random number generator seed

# Output files

Two files are generated the course of a simulation:

1. An output file with simulation results arranged in 10 columns:

- column 1: MC cycle number

- column 2: average membrane separation

- column 3: membrane roughness

- column 4: number of receptor-ligand complexes

- column 5: number of free protein particles

- column 6: total energy of membrane bending

- column 7: total energy of contacts between rafts

- column 8: total energy of particle-raft association

- column 9: total energy of receptor-ligand interactions

- column 10: sum of total energies in columns 6, 7, 8 and 9

2. A configuration file with data arranged in 8 columns:

- column 1: x coordinate

- column 2: y coordinate

- column 3: local position of the upper membrane at this lattice site

- column 4: local position of the lower membrane at this lattice site

- column 5: receptor index (0 if none) at this lattice site

- column 6: ligand index (0 if none) at this lattice site

- column 7: raft index (0 if none) at this lattice site

- column 8: raft index (0 if none) at this lattice site
