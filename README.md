# Monte Carlo code to simulate membrane adhesion

This C++ code can be used to perform Monte Carlo simulations of the lattice-based mesoscale model introduced in:
- Long Li, Jinglei Hu, Bartosz Rozycki, Fan Song, Nano Letters 20(1):722-728 (2020)
- Long Li, Jinglei Hu, Xinghua Shi, Bartosz Rozycki, Fan Song, Soft Matter 17(7):1912-1920 (2021)

Authors: Lukasz Milewski & Bartosz Rozycki
Institute of Physics
Polish Academy of Sciences
October 2023

The development of this code was supported by the National Science Centre via grant number 2021/40/Q/NZ1/00017.

# Usage

Compile the C++ code 

c++ MC_membrane_adhesion.cpp -O -o sim.o

and execute it with 11 parameters at input as in the following example:

./sim.o 0.02 0.1 10.0 0.0 3.0 5.1 10000000 100000000 1000 20000000 1

The input parameters needs to be given in the following order:

1: area concentration of the protein particles (in units of 1/a^2)

2: membrane area fraction of lipid rafts
