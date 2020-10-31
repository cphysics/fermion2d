# QCD fermions in 2d
This repo provides Fortran codes for simulation of 2D gauge theory with fermion content. Dynamics of Fermion in 2 dimension is described by action including dirac operator.

### To run experiment

```
#!/bin/bash
  
  gfortran\
  FERMION.F95\
  Mstart.f95\
  Algebra.f95\
  Printer.f95\
  hdm.f95\
  Calculator.f95\
  Functionary.f95\
  Gamma.f95\
  LFrunner.f95\
  Partition.f95\
  Updator.f95\
  Runner.f95\
  -llapack -lblas

```

### Reference: [Two Dimensional Lattice Gauge Theory with and without Fermion Content](https://digitalcommons.fiu.edu/etd/3224/)
