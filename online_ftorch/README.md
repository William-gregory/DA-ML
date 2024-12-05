## FTorch implementation of an ML-based correction to category sea ice concentrations in SIS2

This repository contains the Python script `NNetwork` which is called in Fortran. The repositories which contain the Fortran code to execute this are in [SIS2](https://github.com/William-gregory/SIS2/tree/ftorch_SPEAR/src).

Specifically, see source code:

- [SIS_slow_thermo.F90](https://github.com/William-gregory/SIS2/tree/ftorch_SPEAR/src/SIS_slow_thermo.F90)
  
  The update to category sea ice concentrations happens after the thermodynamic update to sea ice
- [SIS_ML.F90](https://github.com/William-gregory/SIS2/tree/ftorch_SPEAR/src/SIS_ML.F90)
