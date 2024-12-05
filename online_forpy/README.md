## Forpy implementation of an ML-based correction to category sea ice concentrations in SIS2

This repository contains the Python script `NNetwork` which is called in Fortran. The repositories which contain the Fortran code to execute this are in [SIS2](https://github.com/William-gregory/SIS2/tree/forpy_SPEAR/src).

Specifically, see source code:

- [SIS_slow_thermo.F90](https://github.com/William-gregory/SIS2/tree/forpy_SPEAR/src/SIS_slow_thermo.F90)
  
  The update to category sea ice concentrations happens after the thermodynamic update to sea ice
- [SIS_ML.F90](https://github.com/William-gregory/SIS2/tree/forpy_SPEAR/src/SIS_ML.F90)
  
  State variables are collated into an input array for the network. Part_size is then updated based on the predicted correction, and sea ice variables are subequently updated
  depending on whether ice has been added/removed to a given category
- [Forpy_interface.F90](https://github.com/William-gregory/SIS2/tree/forpy_SPEAR/src/Forpy_interface.F90)
  
  Forpy routine which calls the function inside `NNetwork`
- [forpy_mod.F90](https://github.com/William-gregory/SIS2/tree/forpy_SPEAR/src/forpy_mod.F90)
  
  Forpy interface (unchanged from default version)

For compiling the MOM6/SIS2 model code in order to run Forpy, see the steps outlined in [MOM6-examples](https://github.com/William-gregory/MOM6-examples/blob/forpy_dev/COMPILE_MOM6SIS2.sh).
